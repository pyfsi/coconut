import re
import numpy as np
import os
from collections import OrderedDict

float_pattern = r'[+-]?\d*\.?\d*[eE]?[+-]?\d*'
int_pattern = r'[+-]?\d*'


def get_float(string, keyword):
    result = re.search(keyword + r'\s+\n*(?P<value>' + float_pattern + r')', string)
    if result is None:
        raise RuntimeError(f'keyword not found: {keyword}')

    else:
        return float(result.group('value'))


def get_int(string, keyword):
    result = re.search(keyword + r'\s+\n*(?P<value>' + int_pattern + r')', string)
    if result is None:
        raise RuntimeError(f'keyword not found: {keyword}')

    else:
        return int(result.group('value'))


def get_dict(string, keyword):
    result = re.search(keyword + r'\s+\n*\{.*?\}', string, flags=re.S)
    if result is None:
        raise RuntimeError(f'keyword not found: {keyword}')
    else:
        return result.group()


def get_vector_float_array(string):
    pattern = re.compile(r'\(\s*' + float_pattern + r'\s*' + float_pattern + r'\s*' + float_pattern + r'\s*\)',
                         flags=re.S)
    data_list = re.findall(pattern, string)
    data = np.empty(shape=(len(data_list), 3))
    pattern = re.compile(r'\((.*)\)')
    for i, elem in enumerate(data_list):
        data[i, :] = np.array(re.search(pattern, elem).group(1).split(), dtype=np.float)

    return data


def get_vector_float_array_from_dict(dict_string):
    # keyword = 'value'+r'\s+\n*'+r'nonuniform List\<vector\>'
    keyword = 'value'
    result = re.search(keyword + r'\s+\n*' + 'nonuniform List<vector>\s+\n*\d*\s+\n*' + r'(.*)', dict_string, re.S)

    if result is None:
        raise RuntimeError(f'keyword not found: {keyword}')
    else:

        return get_vector_float_array(result.group(1))


# TODO: To test this function
def get_scalar_float_array(string):
    pattern = re.compile(r'\(\s*' + float_pattern + r'\s*\)',
                         flags=re.S)
    data_list = re.findall(pattern, string)
    data = np.empty(shape=(len(data_list), 3))
    pattern = re.compile(r'\((.*)\)')
    for i, elem in enumerate(data_list):
        data[i, :] = np.array(re.search(pattern, elem).group(1).split(), dtype=np.float)

    return data


def get_boundary_field(file_name, boundary_name):
    with open(file_name, 'r') as f:
        lines = f.read()

    boundary_field_dict = get_dict(lines, boundary_name)
    boundary_field = (get_vector_float_array_from_dict(boundary_field_dict))
    return boundary_field


def get_boundary_points(case_directory, time_folder, boundary_name):
    if time_folder == '0':
        p_file_name = os.path.join(case_directory, 'constant', 'polyMesh/points')
    else:
        p_file_name = os.path.join(case_directory, time_folder, 'polyMesh/points')

    f_file_name = os.path.join(case_directory, 'constant/polyMesh/faces')
    b_file_name = os.path.join(case_directory, 'constant/polyMesh/boundary')

    with open(b_file_name, 'r') as f:
        b_lines = f.read()
        b_dict = get_dict(b_lines, boundary_name)
    nr_faces = get_int(string=b_dict, keyword='nFaces')
    start_face = get_int(string=b_dict, keyword='startFace')

    with open(f_file_name, 'r') as f:
        faces_lines = f.read()
    faces_block = re.search(r'\(.*\)', faces_lines, re.S).group()
    pattern = re.compile(r'(\d*\s*\(.*?\))', flags=re.S)
    faces = re.findall(pattern, faces_block)[start_face:start_face + nr_faces]

    with open(p_file_name, 'r') as f:
        point_lines = f.read()
    point_coords = get_vector_float_array(point_lines)

    pattern = re.compile(r'\((.*)\)')
    point_ids = []
    for i, face in enumerate(faces):
        face_point_ids = re.search(pattern, face).group(1).split()
        for point_id in face_point_ids:
            point_ids.append(int(point_id))
    ordered_points = list(OrderedDict.fromkeys(point_ids))
    return point_coords[ordered_points]


if __name__ == '__main__':
    case_directory = '/lusers/temp/navaneeth/Software/Coconut/coconut/tests/solver_wrappers/openfoam/test_41_pipeCyclic'
    time_folder = '0.1'
    boundary_name = 'walls'
    get_boundary_points(case_directory, time_folder, boundary_name)
