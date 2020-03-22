# This file serves no purpose as of now. When ready, modularize main file using this, similar to analysis functions.

data_load_sources_supported = ['RBPDB (computational)', 'ATTRACT (computational)', 'RBPMap (computational',
                               'POSTAR (experimental)', 'User custom data']
data_load_sources_supported_internal = ["rbpdb", "attract", "rbpmap", "experimental", "custom"]

data_load_sources_functions = {'binding': overall_correlation_analysis, 'per_binding': perBindingSiteAnalysis,
                             'ucsc': ucsc_visualize, 'comp_coop': comp_coop_analysis}
