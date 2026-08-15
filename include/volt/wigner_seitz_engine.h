#pragma once

#include <volt/core/volt.h>
#include <volt/core/simulation_cell.h>
#include <volt/core/particle_property.h>

#include <memory>
#include <vector>

namespace Volt{

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
        bool useMinimumImageConvention = true
    );

    void perform();

    std::shared_ptr<Particles::ParticleProperty> occupancy() const{
        return _occupancy;
    }

    std::shared_ptr<Particles::ParticleProperty> siteIdentifier() const{
        return _siteIdentifier;
    }

    const std::vector<std::size_t>& currentToSite() const{
        return _currentToSite;
    }

    std::size_t vacancyCount() const{ return _vacancyCount; }
    std::size_t interstitialCount() const{ return _interstitialCount; }
    std::size_t antisiteCount() const{ return _antisiteCount; }
    std::size_t occupiedCount() const{ return _occupiedCount; }

private:
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

    std::shared_ptr<Particles::ParticleProperty> _occupancy;
    std::shared_ptr<Particles::ParticleProperty> _siteIdentifier;
    std::vector<std::size_t> _currentToSite;

    std::size_t _vacancyCount = 0;
    std::size_t _interstitialCount = 0;
    std::size_t _antisiteCount = 0;
    std::size_t _occupiedCount = 0;
};

}
