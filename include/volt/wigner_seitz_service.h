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

/**
 * Service wrapper around WignerSeitzEngine. Handles frame validation,
 * particle-property construction and Parquet serialisation of the results.
 */
class WignerSeitzService{
public:
    WignerSeitzService();

    void setAffineMapping(WignerSeitzEngine::AffineMappingType mode);
    void setEliminateCellDeformation(bool enabled);
    void setMinimumImageConvention(bool enabled);
    void setPerTypeOccupancies(bool enabled);

    /**
     * Run the analysis. The reference frame is required and must contain a
     * positive number of atoms. The current and reference frames intentionally
     * may have different atom counts. Two Parquet files are emitted when
     * @p outputBase is non-empty:
     *  - {outputBase}_wigner_seitz.parquet  - per-site / summary defect listing.
     *  - {outputBase}_atoms.parquet         - per-role records (a synthetic
     *    frame of reference sites + extra interstitials), bucket=role drives
     *    the AtomisticExporter grouping (Vacancy/Occupied/Interstitial/Antisite).
     */
    json compute(
        const LammpsParser::Frame& frame,
        const LammpsParser::Frame& referenceFrame,
        const std::string& outputBase = ""
    );

private:
    WignerSeitzEngine::AffineMappingType _affineMapping;
    bool _eliminateCellDeformation;
    bool _useMinimumImageConvention;
    bool _perTypeOccupancies;
};

}
