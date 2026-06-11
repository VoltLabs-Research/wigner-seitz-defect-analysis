#include <volt/wigner_seitz_service.h>
#include <volt/wigner_seitz_engine.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/plugin/output_serializer.h>
#include <volt/utilities/parquet_atom_writer.h>
#include <spdlog/spdlog.h>

#include <chrono>
#include <string>
#include <vector>

namespace Volt{

using namespace Volt::Particles;

namespace{

std::shared_ptr<ParticleProperty> maybeCreateTypeProperty(const LammpsParser::Frame& frame){
    if(frame.types.empty() || frame.natoms <= 0){
        return nullptr;
    }
    auto prop = std::make_shared<ParticleProperty>(
        static_cast<std::size_t>(frame.natoms),
        ParticleProperty::ParticleTypeProperty,
        1,
        true
    );
    int* dst = prop->dataInt();
    const std::size_t n = std::min<std::size_t>(frame.types.size(), prop->size());
    for(std::size_t i = 0; i < n; ++i){
        dst[i] = frame.types[i];
    }
    return prop;
}

std::string roleName(int siteOccupancy, bool isAntisite){
    if(siteOccupancy == 0) return "Vacancy";
    if(isAntisite) return "Antisite";
    if(siteOccupancy == 1) return "Occupied";
    return "Interstitial";
}

}

WignerSeitzService::WignerSeitzService()
    : _affineMapping(WignerSeitzEngine::AffineMappingType::Off),
      _eliminateCellDeformation(false),
      _useMinimumImageConvention(true){}

void WignerSeitzService::setAffineMapping(WignerSeitzEngine::AffineMappingType mode){
    _affineMapping = mode;
}

void WignerSeitzService::setEliminateCellDeformation(bool enabled){
    _eliminateCellDeformation = enabled;
}

void WignerSeitzService::setMinimumImageConvention(bool enabled){
    _useMinimumImageConvention = enabled;
}

json WignerSeitzService::compute(
    const LammpsParser::Frame& frame,
    const LammpsParser::Frame& referenceFrame,
    const std::string& outputBase
){
    auto startTime = std::chrono::high_resolution_clock::now();

    if(frame.natoms <= 0){
        return AnalysisResult::failure("Invalid number of atoms in current frame");
    }
    if(referenceFrame.natoms <= 0){
        return AnalysisResult::failure("Invalid or missing reference frame");
    }

    if(!FrameAdapter::validateSimulationCell(frame.simulationCell)){
        return AnalysisResult::failure("Invalid simulation cell in current frame");
    }
    if(!FrameAdapter::validateSimulationCell(referenceFrame.simulationCell)){
        return AnalysisResult::failure("Invalid simulation cell in reference frame");
    }

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions){
        return AnalysisResult::failure("Failed to create position property for current frame");
    }
    auto refPositions = FrameAdapter::createPositionPropertyShared(referenceFrame);
    if(!refPositions){
        return AnalysisResult::failure("Failed to create position property for reference frame");
    }

    auto types = maybeCreateTypeProperty(frame);
    auto refTypes = maybeCreateTypeProperty(referenceFrame);

    spdlog::info(
        "Starting Wigner-Seitz analysis (ref_sites={}, current_atoms={}, affine={}, mic={})...",
        referenceFrame.natoms,
        frame.natoms,
        static_cast<int>(_affineMapping),
        _useMinimumImageConvention ? "true" : "false"
    );

    WignerSeitzEngine engine(
        positions.get(),
        types.get(),
        frame.simulationCell,
        refPositions.get(),
        refTypes.get(),
        referenceFrame.simulationCell,
        _affineMapping,
        _eliminateCellDeformation,
        _useMinimumImageConvention
    );

    try{
        engine.perform();
    }catch(const std::exception& ex){
        return AnalysisResult::failure(std::string("Wigner-Seitz engine failed: ") + ex.what());
    }

    auto occupancyProp = engine.occupancy();
    const auto& currentToSite = engine.currentToSite();

    const std::size_t nRef = static_cast<std::size_t>(referenceFrame.natoms);
    const std::size_t nCurr = static_cast<std::size_t>(frame.natoms);

    // Build the per-site summary payload (written to _wigner_seitz.parquet).
    json result = AnalysisResult::success();
    result["main_listing"] = {
        { "vacancy_count", engine.vacancyCount() },
        { "interstitial_count", engine.interstitialCount() },
        { "antisite_count", engine.antisiteCount() },
        { "occupied_count", engine.occupiedCount() },
        { "total_sites", nRef },
        { "total_current_atoms", nCurr }
    };

    json perSite = json::array();
    if(occupancyProp){
        const int* occ = occupancyProp->constDataInt();
        for(std::size_t i = 0; i < nRef; ++i){
            const Point3 sitePos = (i < referenceFrame.positions.size())
                ? referenceFrame.positions[i]
                : Point3::Origin();
            const int siteId = (i < referenceFrame.ids.size())
                ? referenceFrame.ids[i]
                : static_cast<int>(i);
            json site;
            site["site_index"] = i;
            site["site_id"] = siteId;
            site["pos"] = { sitePos.x(), sitePos.y(), sitePos.z() };
            site["occupancy"] = occ[i];
            perSite.push_back(site);
        }
    }
    result["per-site-properties"] = perSite;

    // Map each current atom to the reference site it occupies. Reused for both
    // antisite detection and the per-role record construction below.
    std::vector<std::vector<std::size_t>> occupantsBySite(nRef);
    for(std::size_t i = 0; i < nCurr; ++i){
        const std::size_t s = currentToSite[i];
        if(s < nRef){
            occupantsBySite[s].push_back(i);
        }
    }

    // Determine antisite sites: a reference site whose occupancy is exactly
    // one and whose sole occupant has a different particle type.
    std::vector<bool> siteIsAntisite(nRef, false);
    if(!referenceFrame.types.empty() && !frame.types.empty()){
        for(std::size_t s = 0; s < nRef; ++s){
            if(occupantsBySite[s].size() == 1){
                const std::size_t occupant = occupantsBySite[s][0];
                const int refType = (s < referenceFrame.types.size()) ? referenceFrame.types[s] : 0;
                const int currType = (occupant < frame.types.size()) ? frame.types[occupant] : 0;
                if(refType != 0 && currType != 0 && refType != currType){
                    siteIsAntisite[s] = true;
                }
            }
        }
    }

    // Build the per-role record set. Each reference site becomes one record;
    // sites with occupancy >= 2 emit additional Interstitial records (one per
    // extra occupant). These records (NOT the current atoms) become the rows of
    // a SYNTHETIC frame so that vacancies - which have no current atom - are
    // still represented in the _atoms.parquet table. bucket = role drives the
    // four AtomisticExporter structures.
    std::vector<Point3> recPos;
    std::vector<int> recId;
    std::vector<std::string> recRole;
    std::vector<int> recOccupancy;

    // Pass 1: per-site role assignment (Vacancy / Occupied / Antisite / Interstitial).
    if(occupancyProp){
        const int* occ = occupancyProp->constDataInt();
        for(std::size_t s = 0; s < nRef; ++s){
            const int total = occ[s];

            const bool antisiteHere = siteIsAntisite[s];
            const std::string role = roleName(total, antisiteHere);

            // Site position is the reference site, except when a single current
            // atom sits there, in which case we expose its current position
            // (useful for vis tooling following the occupant).
            Point3 pos;
            int atomId;
            if(total == 1 && !occupantsBySite[s].empty()){
                const std::size_t occupant = occupantsBySite[s][0];
                pos = (occupant < frame.positions.size())
                    ? frame.positions[occupant]
                    : Point3::Origin();
                atomId = (occupant < frame.ids.size())
                    ? frame.ids[occupant]
                    : static_cast<int>(occupant);
            }else{
                pos = (s < referenceFrame.positions.size())
                    ? referenceFrame.positions[s]
                    : Point3::Origin();
                atomId = (s < referenceFrame.ids.size())
                    ? referenceFrame.ids[s]
                    : static_cast<int>(s);
            }

            recPos.push_back(pos);
            recId.push_back(atomId);
            recRole.push_back(role);
            recOccupancy.push_back(total);
        }
    }

    // Pass 2: extra interstitial records for every "extra" occupant beyond the
    // first one at sites with total occupancy >= 2.
    for(std::size_t s = 0; s < nRef; ++s){
        const auto& occupants = occupantsBySite[s];
        if(occupants.size() <= 1) continue;
        for(std::size_t k = 1; k < occupants.size(); ++k){
            const std::size_t occupant = occupants[k];
            const Point3 pos = (occupant < frame.positions.size())
                ? frame.positions[occupant]
                : Point3::Origin();
            const int atomId = (occupant < frame.ids.size())
                ? frame.ids[occupant]
                : static_cast<int>(occupant);
            recPos.push_back(pos);
            recId.push_back(atomId);
            recRole.push_back("Interstitial");
            recOccupancy.push_back(static_cast<int>(occupants.size()));
        }
    }

    AnalysisResult::addTiming(result, startTime);

    if(!outputBase.empty()){
        // Synthetic frame: its atoms are the W-S records (ref sites + extra
        // interstitials), not the current atoms. streamAtomsToParquet iterates
        // frame.natoms and reads only positions/ids, so this is what makes
        // vacancies appear in _atoms.parquet.
        LammpsParser::Frame wsFrame;
        wsFrame.natoms = static_cast<int>(recPos.size());
        wsFrame.positions = recPos;
        wsFrame.ids = recId;
        wsFrame.simulationCell = referenceFrame.simulationCell;

        Plugin::serializePluginOutput(outputBase, wsFrame, result, {
            .summaryFileSuffix = "_wigner_seitz", // <base>_wigner_seitz.parquet from `result`
            .bucketResolver = [&recRole](std::size_t i){ return recRole[i]; }, // bucket == role
            .perAtomColumnWriter = [&](ColumnarAtomWriter& w, std::size_t i){
                w.field("occupancy", recOccupancy[i]);
            }
        }); // also writes <base>_atoms.parquet from wsFrame
    }

    spdlog::info(
        "Wigner-Seitz finished: vacancies={}, interstitials={}, antisites={}, occupied_sites={}.",
        engine.vacancyCount(),
        engine.interstitialCount(),
        engine.antisiteCount(),
        engine.occupiedCount()
    );

    return result;
}

}
