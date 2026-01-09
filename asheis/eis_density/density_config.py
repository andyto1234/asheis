"""Configuration for EIS density diagnostics."""

DENSITY_DIAGNOSTICS = {
    'fe_13': {
        'sav_file': 'density_ratios_fe_13_203_82_202_04_.sav',
        'nom_line': 'fe_13_203.83',
        'denom_line': 'fe_13_202.04',
        'vmin': 8,
        'vmax': 10,
        'description': 'Fe XIII 203_202'
    },
    'fe_12_186': {
        'sav_file': 'density_ratios_fe_12_186_88_195_12_.sav',
        'nom_line': 'fe_12_186.88',
        'denom_line': 'fe_12_195.12',
        'vmin': 8,
        'vmax': 10,
        'description': 'Fe XII 186_195'
    },
    'fe_12': {
        'sav_file': 'density_ratios_fe_12_196_64_195_12_.sav',
        'nom_line': 'fe_12_196.64',
        'denom_line': 'fe_12_195.12',
        'vmin': 8,
        'vmax': 10,
        'description': 'Fe XII 196_195'
    },
    'mg_7': {
        'sav_file': 'density_ratios_mg_7_280_74_278_39_.sav',
        'nom_line': 'mg_7_280.74',
        'denom_line': 'mg_7_278.39',
        'vmin': 8,  # Adjust these based on expected density range
        'vmax': 10,  # Adjust these based on expected density range
        'description': 'Mg VII 280_278'
    }
}