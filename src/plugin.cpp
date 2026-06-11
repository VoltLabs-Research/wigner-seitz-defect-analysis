#include <volt/plugin/plugin_main.h>
#include <volt/wigner_seitz_service.h>
#include <volt/wigner_seitz_engine.h>
#include <volt/core/analysis_result.h>
#include <volt/core/lammps_parser.h>

#include <string>

using namespace Volt;

// Wigner-Seitz needs a reference configuration whose site count generally
// DIFFERS from the current atom count (that mismatch is precisely how
// vacancies and interstitials are detected). The framework's
// needsReferenceFrame path in plugin_main.h fatally rejects
// refFrame.natoms != frame.natoms, so we deliberately keep
// needsReferenceFrame = false and parse --reference ourselves below.
static WignerSeitzEngine::AffineMappingType parseAffineMapping(const std::string& s){
    if(s == "toReference" || s == "toReferenceCell") return WignerSeitzEngine::AffineMappingType::ToReference;
    if(s == "toCurrent" || s == "toCurrentCell") return WignerSeitzEngine::AffineMappingType::ToCurrent;
    return WignerSeitzEngine::AffineMappingType::Off; // off | none | noMapping
}

static const Volt::Plugin::PluginDescriptor descriptor{
    .name = "volt-wigner-seitz",
    .description = "Wigner-Seitz Defect Analysis",
    .options = {
        {"--reference", "string", "Reference LAMMPS dump defining the lattice sites (required)", ""},
        {"--affineMapping", "string", "off|toReference|toCurrent", "off"},
        {"--eliminateCellDeformation", "bool", "Eliminate cell deformation before site assignment", "false"},
        {"--minimumImageConvention", "bool", "Use minimum image convention for site assignment", "true"},
        {"--perTypeOccupancies", "bool", "Emit occupancy counts split by particle type", "false"},
    },
    // CHOICE 1: see comment above. The reference frame is parsed manually in
    // the run lambda so the natoms-match guard never runs.
    .needsReferenceFrame = false
};

VOLT_PLUGIN_MAIN(descriptor,
    [](const auto& opts, const LammpsParser::Frame& frame,
       const LammpsParser::Frame* /*refFramePtr: unused, guard bypassed*/,
       const std::string& outputBase) -> Volt::Plugin::json {
        // Resolve the reference dump path. Prefer the modern --reference flag;
        // fall back to the legacy --referenceFrame name for back-compat.
        std::string refFile = CLI::getString(opts, "--reference");
        if(refFile.empty()) refFile = CLI::getString(opts, "--referenceFrame");
        if(refFile.empty()){
            return AnalysisResult::failure(
                "--reference is required for Wigner-Seitz defect analysis");
        }

        LammpsParser::Frame referenceFrame;
        LammpsParser refParser;
        spdlog::info("Parsing reference frame: {}", refFile);
        if(!refParser.parseFile(refFile, referenceFrame)){
            return AnalysisResult::failure("Failed to parse reference frame: " + refFile);
        }
        spdlog::info("Reference loaded: {} sites.", referenceFrame.natoms);

        WignerSeitzService svc;
        svc.setAffineMapping(parseAffineMapping(CLI::getString(opts, "--affineMapping", "off")));
        svc.setEliminateCellDeformation(CLI::getBool(opts, "--eliminateCellDeformation", false));
        svc.setMinimumImageConvention(CLI::getBool(opts, "--minimumImageConvention", true));
        svc.setPerTypeOccupancies(CLI::getBool(opts, "--perTypeOccupancies", false));

        return svc.compute(frame, referenceFrame, outputBase);
    })
