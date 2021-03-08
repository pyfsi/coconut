import re
import numpy as np
import os
from collections import OrderedDict

float_pattern = r'[+-]?\d*\.?\d*[eE]?[+-]?\d*'
int_pattern = r'[+-]?\d+'
delimter = r'[\s\n]+'


def get_float(input_string, keyword):
    result = re.search(keyword + delimter + r'(?P<value>' + float_pattern + r')', input_string)
    if result is None:
        raise RuntimeError(f'keyword not found: {keyword}')

    else:
        return float(result.group('value'))


def get_int(input_string, keyword):
    result = re.search(keyword + delimter + r'(?P<value>' + int_pattern + r')', input_string)
    if result is None:
        raise RuntimeError(f'keyword not found: {keyword}')

    else:
        return int(result.group('value'))


def get_string(input_string, keyword):
    result = re.search(keyword + delimter + r'(?P<value>' + r'\w+' + r')', input_string)
    if result is None:
        raise RuntimeError(f'keyword not found: {keyword}')

    else:
        return result.group('value')


def get_dict(input_string, keyword):
    result = re.search(keyword + delimter + r'\{.*?\}', input_string, flags=re.S)
    if result is None:
        raise RuntimeError(f'keyword not found: {keyword}')
    else:
        return result.group()


def get_vector_array(input_string, is_int=False):
    if is_int:
        pattern = re.compile(
            r'\(' + r'[\s\n]*' + int_pattern + delimter + int_pattern + delimter + int_pattern + r'[\s\n]*\)',
            flags=re.S)
    else:
        pattern = re.compile(
            r'\(' + r'[\s\n]*' + float_pattern + delimter + float_pattern + delimter + float_pattern + r'[\s\n]*\)',
            flags=re.S)
    data_list = re.findall(pattern, input_string)
    data = np.empty(shape=(len(data_list), 3))
    pattern = re.compile(r'\((.*)\)')

    if is_int:
        for i, elem in enumerate(data_list):
            data[i, :] = np.array(re.search(pattern, elem).group(1).strip().split(), dtype=np.int)
    else:
        for i, elem in enumerate(data_list):
            data[i, :] = np.array(re.search(pattern, elem).group(1).strip().split(), dtype=np.float)

    return data


def get_vector_array_from_dict(dict_string, size=None, is_int=False):
    # keyword = 'value'+r'\s+\n*'+r'nonuniform List\<vector\>'
    keyword = 'value'
    result_non_uniform = re.search(
        keyword + delimter + 'nonuniform List<vector>' + delimter + r'\d*' + delimter + r'(.*)', dict_string, re.S)
    result_uniform = re.search(keyword + delimter + r'uniform' + delimter + r'\((.*)\)', dict_string)
    if not result_non_uniform is None:
        return get_vector_array(input_string=result_non_uniform.group(1), is_int=is_int)
    elif not (result_uniform is None or size is None):
        value_list = result_uniform.group(1).strip().split()
        if is_int:
            return np.full((size, len(value_list)), value_list, dtype=np.int)
        else:
            return np.full((size, len(value_list)), value_list, dtype=np.float)
    else:
        raise RuntimeError(f'keyword not found: {keyword}')


def update_vector_array_dict(dict_string, vector_array):
    keyword = 'value'
    result_non_uniform = re.search(
        keyword + delimter + 'nonuniform List<vector>' + delimter + r'(\d*' + delimter + r'\(.*\))', dict_string, re.S)
    result_uniform = re.search(keyword + delimter + r'uniform' + delimter + r'(\(.*\))', dict_string)
    size = vector_array.shape[0]
    replace_string = '\n' + str(size) + '\n' + np.array2string(vector_array, separator=' ', threshold=1e15,
                                                               precision=15).replace('[', '(').replace(
        ']', ')')
    if not result_non_uniform is None:
        return_string = dict_string.replace(result_non_uniform.group(1), replace_string)
        return return_string

    elif not (result_uniform is None):
        return_string = dict_string.replace('uniform', 'nonuniform List<vector>')
        return_string = return_string.replace(result_uniform.group(1), replace_string)
        return return_string

    else:
        raise RuntimeError(f'keyword not found: {keyword}')


def get_scalar_array(input_string, is_int=False):
    scalar_string = re.search(r'\((.*)\)', input_string, re.S).group(1)
    data_list = scalar_string.split()
    if is_int:
        return np.array(data_list, dtype=np.int)
    else:
        return np.array(data_list, dtype=np.float)


def get_scalar_array_from_dict(dict_string, size=None, is_int=False):
    keyword = 'value'
    result_non_uniform = re.search(
        keyword + delimter + 'nonuniform List<scalar>' + delimter + r'\d*' + delimter + r'(.*)', dict_string, re.S)
    result_uniform = re.search(keyword + delimter + r'uniform' + delimter + r'\((.*)\)', dict_string)
    if not result_non_uniform is None:
        return get_scalar_array(input_string=result_non_uniform.group(1), is_int=is_int)
    elif not (result_uniform is None or size is None):
        value_list = result_uniform.group(1).strip().split()
        if is_int:
            return np.full((size, len(value_list)), value_list, dtype=np.int)
        else:
            return np.full((size, len(value_list)), value_list, dtype=np.float)
    else:
        raise RuntimeError(f'keyword not found: {keyword}')


def get_boundary_field(file_name, boundary_name, is_scalar, size=None, is_int=False):
    with open(file_name, 'r') as f:
        file_string = f.read()

    boundary_field_dict = get_dict(input_string=file_string, keyword=boundary_name)

    if is_scalar:
        boundary_field = get_scalar_array_from_dict(dict_string=boundary_field_dict, size=size, is_int=is_int)
    else:
        boundary_field = get_vector_array_from_dict(dict_string=boundary_field_dict, size=size, is_int=is_int)
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
    nr_faces = get_int(input_string=b_dict, keyword='nFaces')
    start_face = get_int(input_string=b_dict, keyword='startFace')

    with open(f_file_name, 'r') as f:
        faces_lines = f.read()
    faces_block = re.search(r'\(.*\)', faces_lines, re.S).group()
    pattern = re.compile(r'(\d*\s*\(.*?\))', flags=re.S)
    faces = re.findall(pattern, faces_block)[start_face:start_face + nr_faces]

    with open(p_file_name, 'r') as f:
        point_file_string = f.read()
    point_coords = get_vector_array(input_string=point_file_string, is_int=False)

    pattern = re.compile(r'\((.*)\)')
    point_ids = []
    for i, face in enumerate(faces):
        face_point_ids = re.search(pattern, face).group(1).split()
        for point_id in face_point_ids:
            point_ids.append(int(point_id))
    unique_point_ids = list(OrderedDict.fromkeys(point_ids))
    return np.array(unique_point_ids), point_coords[unique_point_ids]
