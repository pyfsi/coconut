from sys import argv


def convert(neu='test.neu', inp='base.inp', type='CAX8RH'):

    def print_mod(data, mod):
        ind = 0
        output = ''
        while ind < len(data):
            out = data[ind]
            for i in range(1, min(mod, len(data) - ind)):
                out += ',' + data[ind + i]
            output += out + '\n'
            ind += 8
        return output

    with open(neu, 'r') as in_file:
        with open(inp, 'w') as out_file:
            for line in in_file:
                if 'NUMNP' in line:
                    num_points, num_elems, num_groups, num_boundary_sets, num1, num2 = (int(num) for num in in_file.readline().split())
                    break

            # nodes
            data_nodes = []
            out_file.write('*Part, NAME=PART-1\n*NODE\n')
            for line in in_file:
                if 'NODAL COORDINATES' in line:
                    for _ in range(num_points):
                        num, x, y = tuple(in_file.readline().split())
                        out_file.write(f'{num},{float(x):g},{float(y):g}\n')
                        data_nodes.append(num)
                    break

            # elements
            data_elements = []
            out_file.write(f'*ELEMENT, TYPE={type}\n')
            for line in in_file:
                if 'ELEMENTS/CELLS' in line:
                    for _ in range(num_elems):
                        data = (in_file.readline().split() + in_file.readline().split())
                        num, _, _, p0, p4, p1, p5, p2, p6, p3, p7 = tuple(data)
                        out_file.write(f'{num},{p0},{p1},{p2},{p3},{p4},{p5},{p6},{p7}\n')
                        data_elements.append(num)
                    break

            # groups
            data_groups = {}
            for _ in range(num_groups):
                for line in in_file:
                    if 'ELEMENT GROUP' in line:
                        in_file.readline()
                        name = in_file.readline().strip()
                        in_file.readline()
                        data_groups[name] = []
                        line = in_file.readline()
                        while 'ENDOFSECTION' not in line:
                            data_groups[name] += line.split()
                            line = in_file.readline()
                        break

            # boundary conditions
            data_bc = {}
            for _ in range(num_boundary_sets):
                for line in in_file:
                    if 'BOUNDARY CONDITIONS' in line:
                        bc = in_file.readline().split()
                        name = bc[0]
                        num_nodes = int(bc[2])
                        data_bc[name] = []
                        line = in_file.readline()
                        while 'ENDOFSECTION' not in line:
                            data_bc[name] += line.split()
                            line = in_file.readline()
                        break

            out_file.write(f'*NSET, NSET=NALL\n')
            out_file.write(print_mod(data_nodes, 8))
            for name, data in data_bc.items():
                out_file.write(f'*NSET, NSET={name}\n')
                out_file.write(print_mod(data, 8))
            out_file.write(f'*ELSET, ELSET=EALL\n')
            out_file.write(print_mod(data_elements, 8))
            for name, data in data_groups.items():
                out_file.write(f'*ELSET, ELSET={name}\n')
                out_file.write(print_mod(data, 8))
            out_file.write('*End Part\n')

            out_file.write('*Assembly, NAME=Assembly\n*Instance, NAME=PART-1-1, PART=PART-1\n*End Instance\n')
            out_file.write('*NSET, NSET=NALL, INSTANCE=PART-1-1\n')
            out_file.write(print_mod(data_nodes, 8))
            for name, data in data_bc.items():
                out_file.write(f'*NSET, NSET={name}, INSTANCE=PART-1-1\n')
                out_file.write(print_mod(data, 8))
            out_file.write(f'*ELSET, ELSET=EALL, INSTANCE=PART-1-1\n')
            out_file.write(print_mod(data_elements, 8))
            for name, data in data_groups.items():
                out_file.write(f'*ELSET, ELSET={name}, INSTANCE=PART-1-1\n')
                out_file.write(print_mod(data, 8))
            out_file.write('*End Assembly\n')


convert(argv[1], argv[2], argv[3])
