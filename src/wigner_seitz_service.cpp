#include <volt/wigner_seitz_service.h>
#include <volt/wigner_seitz_engine.h>
#include <volt/wigner_seitz_log.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/plugin/output_serializer.h>
#include <volt/utilities/parquet_atom_writer.h>
#include <volt/core/progress.h>

#include <chrono>
#include <string>
#include <vector>

namespace Volt{

using namespace Volt::Particles;

std::shared_ptr<ParticleProperty> WignerSeitzService::maybeCreateTypeProperty(const LammpsParser::Frame& frame){
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

std::string WignerSeitzService::roleName(int siteOccupancy, bool isAntisite){
    if(siteOccupancy == 0) return "Vacancy";
    if(isAntisite) return "Antisite";
    if(siteOccupancy == 1) return "Occupied";
    return "Interstitial";
}

WignerSeitzService::WignerSeitzService()
    : _affineMapping(AffineMapping::NoMapping),
      _eliminateCellDeformation(false),
      _useMinimumImageConvention(true){}

void WignerSeitzService::setAffineMapping(AffineMapping mode){
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
    WignerSeitzLog& log,
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

    ProgressSteps steps({
        {"Assigning atoms to sites", 1.0},
        {"Building per-site output", 2.0},
        {"Writing outputs", 3.0},
    });
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
        engine.perform(log);
    }catch(const std::exception& ex){
        return AnalysisResult::failure(std::string("Wigner-Seitz engine failed: ") + ex.what());
    }
    log.dump(engine);
    steps.next();

    auto occupancyProp = engine.occupancy();
    const auto& currentToSite = engine.currentToSite();

    const std::size_t nRef = static_cast<std::size_t>(referenceFrame.natoms);
    const std::size_t nCurr = static_cast<std::size_t>(frame.natoms);

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

    std::vector<std::vector<std::size_t>> occupantsBySite(nRef);
    for(std::size_t i = 0; i < nCurr; ++i){
        const std::size_t s = currentToSite[i];
        if(s < nRef){
            occupantsBySite[s].push_back(i);
        }
    }

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

    SiteRecords records;

    const auto pushSite = [&](std::size_t s){
        records.siteIndex.push_back(static_cast<int>(s));
        records.siteId.push_back(s < referenceFrame.ids.size()
            ? referenceFrame.ids[s] : static_cast<int>(s));
        records.siteType.push_back(s < referenceFrame.types.size() ? referenceFrame.types[s] : 0);
    };

    if(occupancyProp){
        const int* occ = occupancyProp->constDataInt();
        for(std::size_t s = 0; s < nRef; ++s){
            const int total = occ[s];

            const bool antisiteHere = siteIsAntisite[s];
            const std::string role = roleName(total, antisiteHere);

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

            records.positions.push_back(pos);
            records.ids.push_back(atomId);
            records.roles.push_back(role);
            records.occupancy.push_back(total);
            pushSite(s);
        }
    }

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
            records.positions.push_back(pos);
            records.ids.push_back(atomId);
            records.roles.push_back("Interstitial");
            records.occupancy.push_back(static_cast<int>(occupants.size()));
            pushSite(s);
        }
    }
    log.dump(records);
    steps.next();

    AnalysisResult::addTiming(result, startTime);

    if(!outputBase.empty()){
        LammpsParser::Frame wsFrame;
        wsFrame.natoms = static_cast<int>(records.positions.size());
        wsFrame.positions = records.positions;
        wsFrame.ids = records.ids;
        wsFrame.simulationCell = referenceFrame.simulationCell;

        Plugin::serializePluginOutput(outputBase, wsFrame, result, {
            .summaryFileSuffix = "_wigner_seitz",
            .bucketResolver = [&records](std::size_t i){ return records.roles[i]; },
            .perAtomColumnWriter = [&records](ColumnarAtomWriter& w, std::size_t i){
                w.field("occupancy", records.occupancy[i]);
                w.field("site_index", records.siteIndex[i]);
                w.field("site_identifier", records.siteId[i]);
                w.field("site_type", records.siteType[i]);
            }
        });
    }

    return result;
}

}
