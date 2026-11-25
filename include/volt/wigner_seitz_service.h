#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>
#include <volt/wigner_seitz_engine.h>

#include <nlohmann/json.hpp>
#include <memory>
#include <string>
#include <vector>

namespace Volt{
using json = nlohmann::json;

class WignerSeitzLog;

class WignerSeitzService{
public:
    // The rows of the atoms table: one per reference site, then one per extra occupant of a site.
    struct SiteRecords{
        std::vector<Point3> positions;
        std::vector<int> ids;
        std::vector<std::string> roles;
        std::vector<int> occupancy;
        std::vector<int> siteIndex;
        std::vector<int> siteId;
        std::vector<int> siteType;
    };

    WignerSeitzService();

    void setAffineMapping(AffineMapping mode);
    void setEliminateCellDeformation(bool enabled);
    void setMinimumImageConvention(bool enabled);

    json compute(
        const LammpsParser::Frame& frame,
        const LammpsParser::Frame& referenceFrame,
        WignerSeitzLog& log,
        const std::string& outputBase = ""
    );

private:
    // The particle types of a frame, or none when it has no type column.
    static std::shared_ptr<Particles::ParticleProperty> maybeCreateTypeProperty(const LammpsParser::Frame& frame);
    // The role a row of the atoms table plays: the site's occupancy, or that its occupant is of another type.
    static std::string roleName(int siteOccupancy, bool isAntisite);

    AffineMapping _affineMapping;
    bool _eliminateCellDeformation;
    bool _useMinimumImageConvention;
};

}
