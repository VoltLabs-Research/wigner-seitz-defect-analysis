#include <volt/plugin/plugin_main.h>
#include <volt/wigner_seitz_service.h>
#include <volt/wigner_seitz_engine.h>
#include <volt/core/analysis_result.h>
#include <volt/core/lammps_parser.h>

#include <string>

using namespace Volt;

static WignerSeitzEngine::AffineMappingType parseAffineMapping(const std::string& s){
    if(s == "toReference" || s == "toReferenceCell") return WignerSeitzEngine::AffineMappingType::ToReference;
    if(s == "toCurrent" || s == "toCurrentCell") return WignerSeitzEngine::AffineMappingType::ToCurrent;
    return WignerSeitzEngine::AffineMappingType::Off;
}

static const Volt::Plugin::PluginDescriptor descriptor{
    .name = "volt-wigner-seitz",
    .description = "Wigner-Seitz Defect Analysis",
    .options = {
        {"--reference", "string", "Reference LAMMPS dump defining the lattice sites (required)", ""},
        {"--affineMapping", "string", "off|toReference|toCurrent", "off"},
        {"--eliminateCellDeformation", "bool", "Eliminate cell deformation before site assignment", "false"},
        {"--minimumImageConvention", "bool", "Use minimum image convention for site assignment", "true"},
    },
    .needsReferenceFrame = false
};

VOLT_PLUGIN_MAIN(descriptor,
    [](const auto& opts, const LammpsParser::Frame& frame,
       const LammpsParser::Frame* /*refFramePtr: unused, guard bypassed*/,
       const std::string& outputBase) -> Volt::Plugin::json {
        std::string refFile = CLI::getString(opts, "--reference");
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

        return svc.compute(frame, referenceFrame, outputBase);
    })
