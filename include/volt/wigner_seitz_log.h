#pragma once

#include <volt/core/run_log.h>
#include <volt/wigner_seitz_service.h>

namespace Volt{

class WignerSeitzEngine;

// The execution log of one Wigner-Seitz analysis. Each phase hands it what it built, and the log writes what that
// holds and how long the phase took. Scripts parse these lines, so their wording stays as it is.
class WignerSeitzLog : public RunLog{
public:
    using RunLog::dump;

    // The site assignment and the defects it found.
    void dump(const WignerSeitzEngine& engine);
    void dump(const WignerSeitzService::SiteRecords& records);
};

}
