from pathlib import Path
import astropy.units as u
def int_to_roman(num):
    """Convert integer to Roman numeral."""
    val = [1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1]
    syb = ["M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I"]
    roman_num = ''
    i = 0
    while num > 0:
        for _ in range(num // val[i]):
            roman_num += syb[i]
            num -= val[i]
        i += 1
    return roman_num

def parse_spectral_line(line):
    """Convert spectral line format (e.g. 'fe_12_195.12') to element and ion key."""
    element, ion_num, _ = line.split(' ')
    element = element.capitalize()
    ion = ion_num
    return (element, ion)

def read_ionization_data(file_path):
    data = {}
    
    with open(file_path, 'r') as file:
        # Skip header lines
        for line in file:
            if not line.startswith('#'):
                break
        
        # Read data lines
        for line in file:
            parts = line.split()
            if len(parts) >= 4:
                element = parts[0]
                ionization = parts[1]
                t_max = float(parts[2])
                v_therm = float(parts[3])
                
                key = (element, ionization)
                data[key] = {'T_MAX': t_max, 'V_THERM': v_therm}

    return data

def get_line_info(line):
    """Get T_MAX and V_THERM for a given spectral line.
    
    Returns:
        dict: Dictionary containing T_MAX and V_THERM values
    """

    file_path = Path(__file__).parent / 'data/eis_width2velocity.dat'
    ionization_data = read_ionization_data(file_path)
    
    key = parse_spectral_line(line)
    
    if key in ionization_data:
        return {
            'T_MAX': ionization_data[key]['T_MAX'],
            'V_THERM': ionization_data[key]['V_THERM'],
            'IONISATION': key
        }
    else:
        raise ValueError(f"No data found for {line}")
    
    import astropy.units as u

def eis_element2mass(line, kg=False):
    element, ionization = parse_spectral_line(line)
    amu2g = 1.6605E-24
    
    data = [
        {"element": "H", "z": 1, "amu": 1.0080},
        {"element": "He", "z": 2, "amu": 4.0026},
        {"element": "Li", "z": 3, "amu": 6.941},
        {"element": "Be", "z": 4, "amu": 9.0122},
        {"element": "B", "z": 5, "amu": 10.811},
        {"element": "C", "z": 6, "amu": 12.0111},
        {"element": "N", "z": 7, "amu": 14.0067},
        {"element": "O", "z": 8, "amu": 15.9994},
        {"element": "F", "z": 9, "amu": 18.9984},
        {"element": "Ne", "z": 10, "amu": 20.179},
        {"element": "Na", "z": 11, "amu": 22.9898},
        {"element": "Mg", "z": 12, "amu": 24.305},
        {"element": "Al", "z": 13, "amu": 26.9815},
        {"element": "Si", "z": 14, "amu": 28.086},
        {"element": "P", "z": 15, "amu": 30.9738},
        {"element": "S", "z": 16, "amu": 32.06},
        {"element": "Cl", "z": 17, "amu": 35.453},
        {"element": "Ar", "z": 18, "amu": 39.948},
        {"element": "K", "z": 19, "amu": 39.102},
        {"element": "Ca", "z": 20, "amu": 40.08},
        {"element": "Sc", "z": 21, "amu": 44.956},
        {"element": "Ti", "z": 22, "amu": 47.90},
        {"element": "V", "z": 23, "amu": 50.9414},
        {"element": "Cr", "z": 24, "amu": 51.996},
        {"element": "Mn", "z": 25, "amu": 54.9380},
        {"element": "Fe", "z": 26, "amu": 55.847},
        {"element": "Co", "z": 27, "amu": 58.9332},
        {"element": "Ni", "z": 28, "amu": 58.71},
        {"element": "Cu", "z": 29, "amu": 63.546},
        {"element": "Zn", "z": 30, "amu": 65.37}
    ]
    
    element = element.upper()
    for item in data:
        if item["element"].upper() == element:
            mass = item["amu"] * amu2g * u.g
            return mass
    
    return -1

