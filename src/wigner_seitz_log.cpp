#include <volt/wigner_seitz_log.h>

#include <volt/wigner_seitz_engine.h>

#include <spdlog/spdlog.h>

#include <cstddef>

namespace Volt{

void WignerSeitzLog::dump(const WignerSeitzEngine& engine){
    const auto occupancy = engine.occupancy();
    const std::size_t sites = occupancy ? occupancy->size() : 0;
    if(occupancy){
        const std::size_t atoms = engine.currentToSite().size();
        const std::size_t assigned = engine.occupiedCount() + engine.interstitialCount();
        spdlog::info("Assigned {} atoms in {}", assigned, phaseElapsed());
        if(assigned < atoms){
            spdlog::warn("{} atoms could not be assigned to a reference site", atoms - assigned);
        }
    }
    spdlog::info(
        "Defects: {} vacancies ({:.3g}% of sites), {} interstitials, {} antisites; {} of {} sites occupied",
        engine.vacancyCount(),
        percentOf(engine.vacancyCount(), sites),
        engine.interstitialCount(),
        engine.antisiteCount(),
        engine.occupiedCount(),
        sites
    );
}

void WignerSeitzLog::dump(const WignerSeitzService::SiteRecords& records){
    spdlog::info("Per-site output ready in {} ({} records)", phaseElapsed(), records.positions.size());
}

}
