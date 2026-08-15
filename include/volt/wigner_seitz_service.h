#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>
#include <volt/wigner_seitz_engine.h>

#include <nlohmann/json.hpp>
#include <memory>
#include <string>

namespace Volt{
using json = nlohmann::json;

class WignerSeitzService{
public:
    WignerSeitzService();

    void setAffineMapping(WignerSeitzEngine::AffineMappingType mode);
    void setEliminateCellDeformation(bool enabled);
    void setMinimumImageConvention(bool enabled);

    json compute(
        const LammpsParser::Frame& frame,
        const LammpsParser::Frame& referenceFrame,
        const std::string& outputBase = ""
    );

private:
    WignerSeitzEngine::AffineMappingType _affineMapping;
    bool _eliminateCellDeformation;
    bool _useMinimumImageConvention;
};

}
