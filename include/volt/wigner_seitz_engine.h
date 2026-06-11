#pragma once

#include <volt/core/volt.h>
#include <volt/core/simulation_cell.h>
#include <volt/core/particle_property.h>

#include <memory>
#include <vector>

namespace Volt{

/**
 * Wigner-Seitz Defect Analysis engine.
 *
 * Counts how many "current" atoms fall inside each Wigner-Seitz cell defined
 * by a reference configuration. Sites whose occupancy is zero are vacancies,
 * sites whose occupancy is greater than one indicate interstitials. When
 * particle-type information is available for both configurations, atoms that
 * fall in a reference site whose type does not match their own are flagged
 * as antisites.
 *
 * The implementation follows the approach used by OVITO's
 * WignerSeitzAnalysisModifier: build a spatial index (k=1 nearest neighbour
 * finder) on the reference positions and assign each current atom to the
 * closest reference site.
 */
class WignerSeitzEngine{
public:
    enum class AffineMappingType{
        Off = 0,
        ToReference,
        ToCurrent
    };

    WignerSeitzEngine(
        Particles::ParticleProperty* positions,
        Particles::ParticleProperty* types,
        const Particles::SimulationCell& cell,
        Particles::ParticleProperty* refPositions,
        Particles::ParticleProperty* refTypes,
        const Particles::SimulationCell& refCell,
        AffineMappingType affineMapping = AffineMappingType::Off,
        bool eliminateCellDeformation = false,
        bool useMinimumImageConvention = true,
        bool perTypeOccupancies = false
    );

    void perform();

    // Occupancy per reference site. componentCount == 1 by default, or
    // ptypeCount when per-type occupancies are requested.
    std::shared_ptr<Particles::ParticleProperty> occupancy() const{
        return _occupancy;
    }

    // Parallel to the reference positions: the site identifier is the index
    // of the reference atom so downstream consumers can re-associate the
    // result with the reference configuration.
    std::shared_ptr<Particles::ParticleProperty> siteIdentifier() const{
        return _siteIdentifier;
    }

    // For each current atom, the index of the reference site it was
    // assigned to (std::numeric_limits<std::size_t>::max() on failure).
    const std::vector<std::size_t>& currentToSite() const{
        return _currentToSite;
    }

    std::size_t vacancyCount() const{ return _vacancyCount; }
    std::size_t interstitialCount() const{ return _interstitialCount; }
    std::size_t antisiteCount() const{ return _antisiteCount; }
    std::size_t occupiedCount() const{ return _occupiedCount; }
    std::size_t ptypeCount() const{ return _ptypeCount; }
    bool perTypeOccupancies() const{ return _perTypeOccupancies; }

private:
    // Build the list of current positions transformed into the reference
    // cell (or into the current cell depending on the affine mapping mode).
    void prepareMappedPositions(std::vector<Point3>& mappedPositions) const;

    Particles::ParticleProperty* _positions = nullptr;
    Particles::ParticleProperty* _types = nullptr;
    Particles::ParticleProperty* _refPositions = nullptr;
    Particles::ParticleProperty* _refTypes = nullptr;

    Particles::SimulationCell _simCell;
    Particles::SimulationCell _simCellRef;

    AffineMappingType _affineMapping = AffineMappingType::Off;
    bool _eliminateCellDeformation = false;
    bool _useMinimumImageConvention = true;
    bool _perTypeOccupancies = false;

    std::shared_ptr<Particles::ParticleProperty> _occupancy;
    std::shared_ptr<Particles::ParticleProperty> _siteIdentifier;
    std::vector<std::size_t> _currentToSite;

    std::size_t _vacancyCount = 0;
    std::size_t _interstitialCount = 0;
    std::size_t _antisiteCount = 0;
    std::size_t _occupiedCount = 0;
    std::size_t _ptypeCount = 0;
};

}
