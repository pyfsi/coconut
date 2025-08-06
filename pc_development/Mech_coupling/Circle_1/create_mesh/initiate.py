import math
import re

# Replace each placeholder |name| with its value
def replace_placeholders(text, values):
    def repl(match):
        key = match.group(1)
        return str(values.get(key, f'<undefined:{key}>'))
    return re.sub(r'\|(\w+)\|', repl, text)

# Define parameters and any expressions
params = {
    'L': 0.1,
    'H': 0.1,
    'R': 0.01,
    'Xc': 0.05,
    'Yc': 0.05,
    'Xintervals': 50,
    'Yintervals': 50,
    'ratio': 1.064,
    'DiagonalIntervals': 60,
}

# Trig expressions for arc vertices at 45°, 135°, etc.
pi = math.pi
params['vx1'] = params['Xc'] + params['R'] * math.cos(pi / 4)
params['vy1'] = params['Yc'] + params['R'] * math.sin(pi / 4)

params['vx2'] = params['Xc'] + params['R'] * math.cos(3 * pi / 4)
params['vy2'] = params['Yc'] + params['R'] * math.sin(3 * pi / 4)

params['vx3'] = params['Xc'] + params['R'] * math.cos(5 * pi / 4)
params['vy3'] = params['Yc'] + params['R'] * math.sin(5 * pi / 4)

params['vx4'] = params['Xc'] + params['R'] * math.cos(7 * pi / 4)
params['vy4'] = params['Yc'] + params['R'] * math.sin(7 * pi / 4)

# Read the template script with placeholders
with open('mesh.jou', 'r') as f:
    template = f.read()

# Processed output
final_script = replace_placeholders(template, params)

# Write to output file
output_filename = 'mesh_filled.jou'
with open(output_filename, 'w') as f:
    f.write(final_script)

# Final check for unresolved placeholders
undefined_placeholders = re.findall(r'<undefined:([^>]+)>', final_script)
if undefined_placeholders:
    raise ValueError(f'The following parameters were not defined: {set(undefined_placeholders)}')

print(f'GAMBIT script written to: {output_filename}')
