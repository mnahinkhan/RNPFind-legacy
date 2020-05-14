from colors import orange, blue, dark_green
from rbpdb_data_load import rbpdb_data_load, rbpdb_column_names, rbpdb_default_label_index, rbpdb_columns_of_interest, \
    rbpdb_default_mouse_over_index, rbpdb_column_descriptions
from attract_data_load import attract_data_load, attract_column_names, attract_default_label_index, \
    attract_default_mouse_over_index, attract_columns_of_interest, attract_descriptions
from postar_data_load import postar_data_load, postar_column_names, postar_default_label_index, \
    postar_columns_of_interest, postar_default_mouse_over_index, postar_column_descriptions
from custom_data_load import custom_data_load

data_load_sources_supported_short = ["rbpdb", "attract", "rbpmap", "postar", "custom"]
data_load_sources_supported_long = ['RBPDB (computational)', 'ATTRACT (computational)', 'RBPMap (computational',
                                    'POSTAR (experimental)', 'User custom data']

data_load_sources_functions = {'rbpdb': rbpdb_data_load, 'attract': attract_data_load, 'postar': postar_data_load,
                               "custom": custom_data_load}

column_data = {"postar": {"names": postar_column_names, "default_label": postar_default_label_index,
                          "interest": postar_columns_of_interest,
                          "default_mouse_over": postar_default_mouse_over_index,
                          "descriptions": postar_column_descriptions},

               "attract": {"names": attract_column_names, "default_label": attract_default_label_index,
                           "interest": attract_columns_of_interest,
                           "default_mouse_over": attract_default_mouse_over_index,
                           "descriptions": attract_descriptions},

               "rbpdb": {"names": rbpdb_column_names, "default_label": rbpdb_default_label_index,
                         "interest": rbpdb_columns_of_interest,
                         "default_mouse_over": rbpdb_default_mouse_over_index,
                         "descriptions": rbpdb_column_descriptions}
               }

data_load_source_colors = {"postar": orange, "attract": blue, "rbpdb": dark_green}
