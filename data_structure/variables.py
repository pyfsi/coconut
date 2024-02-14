variables_dimensions = {
                        'pressure': 1,
                        'traction': 3,
                        'displacement': 3,
                        'density': 1,
                        'area': 1,
                        'heat_flux': 1,
                        'temperature': 1
                        }

# each variable key has a list as value containing a prefix for the corresponding file and an integer related to the dimensions of the file
accepted_variables = {
    'in': {
        'displacement': ['nodes_update', 3],
        'temperature': ['temperature', 1],
        'heat_flux': ['heat_flux', 0]
    },
    'out':{
        'pressure': ['pressure_traction', 0],
        'traction': ['repeat', 0], # catch variable stored in the same file to avoid double reading of the same file
        'temperature': ['temperature', 1],
        'heat_flux': ['heat_flux', 0]
    }
}
