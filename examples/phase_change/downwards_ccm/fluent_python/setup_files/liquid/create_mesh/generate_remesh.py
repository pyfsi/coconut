import re
import math

# Input and Output filenames
INPUT_FILE = "identified_edges.trn"
OUTPUT_FILE = "remesh.jou"


def parse_vector(line):
    """Parses (x, y, z) from a line containing 'Centroid: ( x, y, z )'"""
    match = re.search(r'\(\s*([0-9.eE+-]+),\s*([0-9.eE+-]+),\s*([0-9.eE+-]+)\)', line)
    if match:
        return [float(match.group(1)), float(match.group(2)), float(match.group(3))]
    return None


def parse_gambit_trn(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    edges = []
    current_edge = {}

    # State flags
    reading_faces = False
    reading_zone = False

    for line in lines:
        line = line.strip()

        # 1. Detect Edge Start
        if line.startswith("Edge:"):
            # Save previous edge if exists
            if current_edge:
                edges.append(current_edge)
            # Start new edge
            edge_name = line.split(":")[1].strip()
            current_edge = {
                'name': edge_name,
                'faces': [],
                'zone': None,
                'first_size': 0.0,
                'last_size': 0.0,
                'centroid': [0, 0, 0],
                'length': 0.0  # Added field for Arc Length
            }
            reading_faces = False
            reading_zone = False
            continue

        if not current_edge:
            continue

        # 2. Parse Sizes (for Reversal Logic)
        if "First/Last Size:" in line:
            parts = line.split(":")[-1].strip().split("/")
            if len(parts) == 2:
                current_edge['first_size'] = float(parts[0])
                current_edge['last_size'] = float(parts[1])

        # 3. Parse Faces (for deletion command)
        if "Edge used by the following" in line:
            reading_faces = True
            continue
        if reading_faces:
            if line.startswith("Edge contained") or line.startswith("Centroid") or line.startswith(
                    "-") or line.startswith("Arc Length"):
                reading_faces = False
            else:
                if line: current_edge['faces'].append(line)

        # 4. Parse Boundary Zone
        if "Edge contained in the following" in line:
            reading_zone = True
            continue
        if reading_zone:
            current_edge['zone'] = line
            reading_zone = False

        # 5. Parse Centroid (for Linking Logic)
        if "Centroid:" in line:
            current_edge['centroid'] = parse_vector(line)

        # 6. Parse Arc Length (NEW)
        if "Arc Length:" in line:
            # Format: Arc Length: 0.000118
            try:
                val = line.split(":")[1].strip()
                current_edge['length'] = float(val)
            except (IndexError, ValueError):
                pass

    # Append the last edge
    if current_edge:
        edges.append(current_edge)

    return edges


def get_distance(p1, p2):
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2)


def generate_journal(edges):
    # --- Categorize edges ---
    internal_edges = [e for e in edges if e['zone'] == 'internal']
    outer_edges = [e for e in edges if e['zone'] == 'outer']
    inner_edges = [e for e in edges if e['zone'] == 'inner']

    # --- Collect all unique faces for deletion ---
    all_faces = set()
    for e in edges:
        for face in e['faces']:
            all_faces.add(face)
    all_faces.add("fluid.com")
    sorted_faces = sorted(list(all_faces))

    # --- Logic 1: Internal Edges (Reversal) ---
    internal_normal = []
    internal_reverse = []

    for e in internal_edges:
        if e['first_size'] > e['last_size']:
            internal_reverse.append(e['name'])
        else:
            internal_normal.append(e['name'])

    # --- Logic 2: Split Inner Edges by Length (NEW) ---
    # Threshold chosen based on your data: Long ~0.05, Short ~0.0006
    # A threshold of 0.01 safely separates them.
    LENGTH_THRESHOLD = 0.005

    inner_long = [e['name'] for e in inner_edges if e['length'] > LENGTH_THRESHOLD]
    inner_short = [e['name'] for e in inner_edges if e['length'] <= LENGTH_THRESHOLD]

    # --- Logic 3: Linking Outer to Inner ---
    links = []
    used_inner = set()

    for out_e in outer_edges:
        best_match = None
        min_dist = float('inf')

        for in_e in inner_edges:
            dist = get_distance(out_e['centroid'], in_e['centroid'])
            if dist < min_dist:
                min_dist = dist
                best_match = in_e['name']

        if best_match:
            links.append((out_e['name'], best_match))

    # --- Write Journal ---
    with open(OUTPUT_FILE, 'w') as f:
        f.write("/LIQUID DOMAIN\n\n")
        f.write('import uns "case_restart.cas" d2\n')

        # Face Delete
        if sorted_faces:
            f.write('face delete ' + ' '.join([f'"{face}"' for face in sorted_faces]) + ' onlymesh lowertopology\n\n')

        # Internal Meshing
        if internal_normal:
            names = ' '.join([f'"{name}"' for name in internal_normal])
            f.write(f'edge mesh {names} successive ratio1 1.2 intervals 9\n')

        if internal_reverse:
            names = ' '.join([f'"{name}"' for name in internal_reverse])
            f.write(f'edge mesh {names} successive ratio1 1.2 intervals 9 reverse\n')

        # --- CHANGED: Inner Meshing (Split by Length) ---
        if inner_long:
            long_names = ' '.join([f'"{name}"' for name in inner_long])
            f.write(f'edge mesh {long_names} size 6.25e-5\n')

        if inner_short:
            short_names = ' '.join([f'"{name}"' for name in inner_short])
            f.write(f'edge mesh {short_names} size 3e-5\n')

        # --- Linking (Link Outer TO Inner) ---
        for src, trg in links:
            f.write(f'edge link "{src}" "{trg}"\n')

        f.write('\n')

        # Face Meshing
        if sorted_faces:
            f.write('face mesh ' + ' '.join([f'"{face}"' for face in sorted_faces]) + ' submap\n')

        # Physics Modify
        if internal_edges:
            all_int_names = ' '.join([f'"{e["name"]}"' for e in internal_edges])
            f.write(f'physics modify "internal" btype "INTERIOR" edge {all_int_names}\n')

        if outer_edges:
            outer_names_str = ' '.join([f'"{e["name"]}"' for e in outer_edges])
            f.write(f'physics modify "outer" btype "WALL" edge {outer_names_str}\n')

        if inner_edges:
            inner_names_str = ' '.join([f'"{e["name"]}"' for e in inner_edges])
            f.write(f'physics modify "inner" btype "WALL" edge {inner_names_str}\n')

        f.write('\nexport fluent5 "remesh.msh" nozval\n')
        f.write('save\n')
        f.write('end\n')


if __name__ == "__main__":
    try:
        data = parse_gambit_trn(INPUT_FILE)
        generate_journal(data)
        print(f"Successfully generated {OUTPUT_FILE}")
    except Exception as e:
        print(f"Error: {e}")