#include <volt/plugin/plugin_main.h>
#include <volt/wigner_seitz_service.h>
#include <volt/wigner_seitz_engine.h>
#include <volt/wigner_seitz_log.h>
#include <volt/core/analysis_result.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/progress.h>

#include <string>

using namespace Volt;

static const Volt::Plugin::PluginDescriptor descriptor{
    .name = "volt-wigner-seitz",
    .description = "Wigner-Seitz Defect Analysis",
    .options = {
        {"--reference", "string", "Reference LAMMPS dump defining the lattice sites (required)", ""},
        {"--affineMapping", "string", "off|toReference|toCurrent", "off"},
        {"--eliminateCellDeformation", "bool", "Eliminate cell deformation before site assignment", "false"},
        {"--minimumImageConvention", "bool", "Use minimum image convention for site assignment", "true"},
    },
    .needsReferenceFrame = false,
    .referenceAtomCount = Volt::Plugin::ReferenceAtomCount::Any
};

VOLT_PLUGIN_MAIN(descriptor,
    [](const auto& opts, const LammpsParser::Frame& frame,
       const LammpsParser::Frame* referenceFrame,
       const std::string& outputBase) -> Volt::Plugin::json {
        if(!referenceFrame){
            return AnalysisResult::failure(
                "--reference is required for Wigner-Seitz defect analysis");
        }

        // pluginMain has read the reference; its step stays so that the job's progress values do not move.
        ProgressSteps steps({
            {"Reading the reference frame", 1.0},
            {"Running the Wigner-Seitz analysis", 5.0},
        });
        steps.next();
        WignerSeitzLog log;

        const bool eliminateCellDeformation = CLI::getBool(opts, "--eliminateCellDeformation", false);
        if(eliminateCellDeformation){
            spdlog::warn("eliminateCellDeformation has no effect in this implementation; affineMapping alone controls cell mapping");
        }

        WignerSeitzService svc;
        svc.setAffineMapping(parseAffineMapping(CLI::getString(opts, "--affineMapping", "off")).value_or(AffineMapping::NoMapping));
        svc.setEliminateCellDeformation(eliminateCellDeformation);
        svc.setMinimumImageConvention(CLI::getBool(opts, "--minimumImageConvention", true));

        return svc.compute(frame, *referenceFrame, log, outputBase);
    })
