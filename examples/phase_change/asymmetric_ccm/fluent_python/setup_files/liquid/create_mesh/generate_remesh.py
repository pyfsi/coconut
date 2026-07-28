import re
import math

# Input and Output filenames
INPUT_FILE = "identified_edges.trn"
OUTPUT_FILE = "remesh.jou"


def parse_vector(line):
    """Parses (x, y, z) from a line containing 'Centroid: ( x, y, z )'"""
    # Regex to find numbers inside parenthesis
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
                'centroid': [0, 0, 0]
            }
            reading_faces = False
            reading_zone = False
            continue

        if not current_edge:
            continue

        # 2. Parse Sizes (for Reversal Logic)
        if "First/Last Size:" in line:
            # Format: First/Last Size: 0.000045/0.000010
            parts = line.split(":")[-1].strip().split("/")
            if len(parts) == 2:
                current_edge['first_size'] = float(parts[0])
                current_edge['last_size'] = float(parts[1])

        # 3. Parse Faces (for deletion command)
        if "Edge used by the following" in line:
            reading_faces = True
            continue
        if reading_faces:
            if line.startswith("Edge contained") or line.startswith("Centroid") or line.startswith("-"):
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

    # Append the last edge
    if current_edge:
        edges.append(current_edge)

    return edges


def get_distance(p1, p2):
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2)


def generate_journal(edges):
    # --- categorize edges ---
    internal_edges = [e for e in edges if e['zone'] == 'internal']
    outer_edges = [e for e in edges if e['zone'] == 'outer']
    inner_edges = [e for e in edges if e['zone'] == 'inner']

    # --- Collect all unique faces for deletion ---
    all_faces = set()
    for e in edges:
        for face in e['faces']:
            all_faces.add(face)
    # Ensure fluid.com is included if not picked up
    all_faces.add("fluid.com")
    # Sort for consistent output
    sorted_faces = sorted(list(all_faces))

    # --- Logic 1: Internal Edges (Reversal) ---
    internal_normal = []
    internal_reverse = []

    for e in internal_edges:
        if e['first_size'] > e['last_size']:
            internal_reverse.append(e['name'])
        else:
            internal_normal.append(e['name'])

    # --- Logic 2: Linking Outer to Inner ---
    # Find closest pairs
    links = []  # List of tuples (source=outer, target=inner)

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
        f.write('face delete ' + ' '.join([f'"{face}"' for face in sorted_faces]) + ' onlymesh lowertopology\n\n')

        # Internal Meshing
        if internal_normal:
            names = ' '.join([f'"{name}"' for name in internal_normal])
            f.write(f'edge mesh {names} successive ratio1 1.16 intervals 14\n')

        if internal_reverse:
            names = ' '.join([f'"{name}"' for name in internal_reverse])
            f.write(f'edge mesh {names} successive ratio1 1.16 intervals 14 reverse\n')

        # --- CHANGED: Inner Meshing (Mesh Inner first) ---
        inner_names = ' '.join([f'"{e["name"]}"' for e in inner_edges])
        f.write(f'edge mesh {inner_names} size 8e-5\n')

        # --- Linking (Link Outer TO Inner) ---
        # The 'links' list contains (Outer, Inner).
        # GAMBIT syntax: edge link "EdgeToChange" "ReferenceEdge"
        # So we link "Outer" (Unmeshed) "Inner" (Meshed)
        for src, trg in links:
            f.write(f'edge link "{src}" "{trg}"\n')

        f.write('\n')

        # Face Meshing
        f.write('face mesh ' + ' '.join([f'"{face}"' for face in sorted_faces]) + ' submap\n')

        # Physics Modify
        all_int_names = ' '.join([f'"{e["name"]}"' for e in internal_edges])
        f.write(f'physics modify "internal" btype "INTERIOR" edge {all_int_names}\n')

        # NOTE: Ensure explicit naming for outer/inner so they aren't lost in next step
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