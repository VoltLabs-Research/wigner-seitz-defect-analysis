#include <volt/wigner_seitz_engine.h>
#include <volt/analysis/nearest_neighbor_finder.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace Volt{

using namespace Particles;

namespace{

constexpr std::size_t kInvalidSite = std::numeric_limits<std::size_t>::max();

// Minimal ParticleProperty wrapper that owns a Point3 buffer. Used so we
// can feed transformed positions to the NearestNeighborFinder, which expects
// a PositionProperty regardless of whether the data belongs to the
// reference or current configuration.
std::shared_ptr<ParticleProperty> makePositionProperty(const std::vector<Point3>& source){
    auto prop = std::make_shared<ParticleProperty>(
        source.size(),
        ParticleProperty::PositionProperty,
        3,
        false
    );
    Point3* dst = prop->dataPoint3();
    if(dst && !source.empty()){
        std::memcpy(dst, source.data(), source.size() * sizeof(Point3));
    }
    return prop;
}

}

WignerSeitzEngine::WignerSeitzEngine(
    ParticleProperty* positions,
    ParticleProperty* types,
    const SimulationCell& cell,
    ParticleProperty* refPositions,
    ParticleProperty* refTypes,
    const SimulationCell& refCell,
    AffineMappingType affineMapping,
    bool eliminateCellDeformation,
    bool useMinimumImageConvention,
    bool perTypeOccupancies
)
    : _positions(positions),
      _types(types),
      _refPositions(refPositions),
      _refTypes(refTypes),
      _simCell(cell),
      _simCellRef(refCell),
      _affineMapping(affineMapping),
      _eliminateCellDeformation(eliminateCellDeformation),
      _useMinimumImageConvention(useMinimumImageConvention),
      _perTypeOccupancies(perTypeOccupancies){}

void WignerSeitzEngine::prepareMappedPositions(std::vector<Point3>& mappedPositions) const{
    const std::size_t nCurr = _positions->size();
    mappedPositions.resize(nCurr);
    const Point3* src = _positions->constDataPoint3();

    // Affine mapping policy mirrors OVITO's WignerSeitzAnalysisModifier.
    //  - Off:          copy positions verbatim.
    //  - ToReference:  reduced in current cell -> absolute in reference cell.
    //  - ToCurrent:    reduced in current cell -> absolute in current cell,
    //                  meaning no transform on current atoms; the reference
    //                  cell is the one that gets mapped when the finder is
    //                  built (handled by the caller when applicable).
    const bool mapIntoReference = (_affineMapping == AffineMappingType::ToReference);

    if(!mapIntoReference){
        for(std::size_t i = 0; i < nCurr; ++i){
            mappedPositions[i] = src[i];
        }
        return;
    }

    // Mapping into the reference cell. When "eliminate cell deformation" is
    // on, the reference cell matrix replaces the current one for the
    // reduced-to-absolute step, which effectively normalises away any
    // uniform cell deformation before assigning sites.
    const AffineTransformation currInv = _simCell.inverseMatrix();
    const AffineTransformation refMat = _simCellRef.matrix();

    for(std::size_t i = 0; i < nCurr; ++i){
        Point3 reduced = currInv * src[i];
        if(_useMinimumImageConvention){
            const auto pbc = _simCellRef.pbcFlags();
            for(std::size_t k = 0; k < 3; ++k){
                if(pbc[k]){
                    reduced[k] -= std::floor(reduced[k]);
                }
            }
        }
        mappedPositions[i] = refMat * reduced;
    }

    (void)_eliminateCellDeformation; // Reserved for future per-site remapping.
}

void WignerSeitzEngine::perform(){
    if(!_positions || !_refPositions){
        throw std::runtime_error("WignerSeitzEngine: null input properties.");
    }

    const std::size_t nCurr = _positions->size();
    const std::size_t nRef = _refPositions->size();

    _vacancyCount = 0;
    _interstitialCount = 0;
    _antisiteCount = 0;
    _occupiedCount = 0;
    _currentToSite.assign(nCurr, kInvalidSite);

    if(nRef == 0){
        _occupancy.reset();
        _siteIdentifier.reset();
        _vacancyCount = 0;
        return;
    }

    // Determine the number of particle-type columns when per-type
    // occupancies are requested.
    std::size_t ptypeCount = 1;
    if(_perTypeOccupancies && _types){
        int maxType = 0;
        const int* t = _types->constDataInt();
        for(std::size_t i = 0; i < nCurr; ++i){
            if(t[i] > maxType) maxType = t[i];
        }
        if(_refTypes){
            const int* rt = _refTypes->constDataInt();
            for(std::size_t i = 0; i < nRef; ++i){
                if(rt[i] > maxType) maxType = rt[i];
            }
        }
        ptypeCount = std::max<std::size_t>(1, static_cast<std::size_t>(maxType));
    }
    _ptypeCount = ptypeCount;

    // Use the explicit (count, DataType, componentCount, stride, init) ctor
    // because UserProperty maps to DataType::Void which rejects components>0.
    _occupancy = std::make_shared<ParticleProperty>(
        nRef, Particles::DataType::Int, ptypeCount, 0, true
    );
    _occupancy->setType(ParticleProperty::UserProperty);
    _siteIdentifier = std::make_shared<ParticleProperty>(
        nRef, Particles::DataType::Int64, 1, 0, true
    );
    _siteIdentifier->setType(ParticleProperty::IdentifierProperty);

    // Fill site identifiers with the reference atom index (0..nRef-1) so
    // downstream consumers keep a stable reference.
    for(std::size_t i = 0; i < nRef; ++i){
        _siteIdentifier->setInt64(i, static_cast<std::int64_t>(i));
    }

    // The NearestNeighborFinder operates on a ParticleProperty + a
    // SimulationCell. Depending on affineMapping we either:
    //   - build it on the reference positions/cell (Off, ToReference), or
    //   - build it on reference positions remapped into the current cell
    //     (ToCurrent).
    std::vector<Point3> mappedRef;
    ParticleProperty* finderPositions = _refPositions;
    SimulationCell finderCell = _simCellRef;
    std::shared_ptr<ParticleProperty> ownedRefPositions;

    if(_affineMapping == AffineMappingType::ToCurrent){
        mappedRef.resize(nRef);
        const Point3* refPtr = _refPositions->constDataPoint3();
        const AffineTransformation refInv = _simCellRef.inverseMatrix();
        const AffineTransformation currMat = _simCell.matrix();
        const auto pbc = _simCell.pbcFlags();
        for(std::size_t i = 0; i < nRef; ++i){
            Point3 reduced = refInv * refPtr[i];
            if(_useMinimumImageConvention){
                for(std::size_t k = 0; k < 3; ++k){
                    if(pbc[k]){
                        reduced[k] -= std::floor(reduced[k]);
                    }
                }
            }
            mappedRef[i] = currMat * reduced;
        }
        ownedRefPositions = makePositionProperty(mappedRef);
        finderPositions = ownedRefPositions.get();
        finderCell = _simCell;
    }

    NearestNeighborFinder finder(1);
    if(!finder.prepare(finderPositions, finderCell)){
        throw std::runtime_error("WignerSeitzEngine: failed to build nearest-neighbor finder on reference sites.");
    }

    // Prepare mapped current positions.
    std::vector<Point3> mappedCurrent;
    prepareMappedPositions(mappedCurrent);

    // Atomic counters driving thread-safe occupancy updates.
    std::vector<std::atomic<int>> atomicOccupancy(nRef * ptypeCount);
    for(auto& c : atomicOccupancy) c.store(0, std::memory_order_relaxed);

    std::atomic<std::size_t> antisiteLocal{0};

    const int* currTypes = (_types ? _types->constDataInt() : nullptr);
    const int* refTypesPtr = (_refTypes ? _refTypes->constDataInt() : nullptr);

    tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, nCurr, 1024),
        [&](const tbb::blocked_range<std::size_t>& r){
            // CoreToolkit instantiates Query<{4,16,18,32,64,128}>; use 4 and
            // take only the nearest result.
            NearestNeighborFinder::Query<4> q(finder);
            for(std::size_t i = r.begin(); i < r.end(); ++i){
                // includeSelf=true is REQUIRED: with includeSelf=false the finder
                // silently skips zero-distance hits, so a current atom sitting
                // exactly on a reference site (the common case when the current
                // config is a perfect lattice minus defects) would be assigned
                // to the second-nearest site instead of its correct site.
                q.findNeighbors(mappedCurrent[i], /*includeSelf=*/true);
                const auto& res = q.results();
                if(res.size() == 0) continue;

                const std::size_t siteIndex = res[0].index;
                _currentToSite[i] = siteIndex;

                std::size_t component = 0;
                if(_perTypeOccupancies && currTypes){
                    const int t = currTypes[i];
                    if(t >= 1 && static_cast<std::size_t>(t) <= ptypeCount){
                        component = static_cast<std::size_t>(t - 1);
                    }
                }
                atomicOccupancy[siteIndex * ptypeCount + component]
                    .fetch_add(1, std::memory_order_relaxed);

                if(refTypesPtr && currTypes){
                    if(refTypesPtr[siteIndex] != currTypes[i]){
                        antisiteLocal.fetch_add(1, std::memory_order_relaxed);
                    }
                }
            }
        }
    );

    // Commit atomic occupancy into the output property.
    int* occData = _occupancy->dataInt();
    for(std::size_t i = 0; i < nRef * ptypeCount; ++i){
        occData[i] = atomicOccupancy[i].load(std::memory_order_relaxed);
    }

    // Tally site statistics.
    for(std::size_t i = 0; i < nRef; ++i){
        int total = 0;
        for(std::size_t c = 0; c < ptypeCount; ++c){
            total += occData[i * ptypeCount + c];
        }
        if(total == 0){
            ++_vacancyCount;
        }else{
            ++_occupiedCount;
            if(total > 1){
                _interstitialCount += static_cast<std::size_t>(total - 1);
            }
        }
    }

    _antisiteCount = antisiteLocal.load(std::memory_order_relaxed);
}

}
