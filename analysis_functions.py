from overAllCorrAnalysis import overall_correlation_analysis
from perBindingSiteAnalysis import perBindingSiteAnalysis
from ucsc_visualize import ucsc_visualize
from competition_cooperation_analysis import comp_coop_analysis

analysis_methods_supported_short = ['binding', 'per_binding', 'ucsc', 'comp_coop']

analysis_methods_supported_long = ["Binding correlation analysis", "Per-binding-site analysis",
                                   "Visualize on UCSC Genome Browser", "Competition-Cooperation Visualization"]

analysis_method_functions = {'binding': overall_correlation_analysis, 'per_binding': perBindingSiteAnalysis,
                             'ucsc': ucsc_visualize, 'comp_coop': comp_coop_analysis}
