import re

coco_abaqus_step_pattern = r'coconut_abaqusstep\n'
abq_keywords = ['step', 'boundary conditions',
            'interactions', 'loads',
            'controls', 'output requests',
            'field output', 'history output',
]
start_pattern = [r'\*\* '+key+r'.*[\s]+' for key in abq_keywords]
abq_pattern_dict = {key: pattern for key, pattern in zip(abq_keywords, start_pattern)}
borderline = r'\*\* -+\n'

def read_input_file(input_file):
    """
    Reads the input file and stores the information in a dictionary.
    Assumes the input file contains keywords commonly found in ABAQUS .inp files.
    """
    abaqus_case_settings = {"header": None,}
    for key in abq_keywords:
        abaqus_case_settings[key] = ""

    # read input file and store in data as string
    with open(input_file, 'r') as file:
        data = file.read()

    # check if coconut_abaqusstep is last step in input file
    coco_step_idx = [n.span() for n in re.finditer(coco_abaqus_step_pattern, data, re.I)]
    if not coco_step_idx:
        coco_step_idx = [n.span() for n in re.finditer(r"step-1", data, re.I)]
    if not coco_step_idx:
        raise RuntimeError("No step definition with the name coco_abaqusstep nor step-1 has been found.")
    borderline_idxs = [n.span() for n in re.finditer(borderline, data,)]
    is_coco_step_lastline = coco_step_idx[-1][-1] > borderline_idxs[-1][-1]
    if not is_coco_step_lastline:
        raise RuntimeError("coconut_abaqusStep should be defined as the last step of the ABAQUS input file.")

    # ignore lines before coconut_abaqusstep. It should be the last step.
    coco_step_data = data[borderline_idxs[-1][0]:]
    # print(coco_step_data)

    # get indices of found keywords
    keyword_idxs = {}
    for key in abq_keywords:
        indexes = [[n.span()[0], n.span()[1]] for n in re.finditer(abq_pattern_dict[key], coco_step_data, re.I)]
        keyword_idxs[key] = indexes

    # extract sorted start and end indices
    idx_entries = {"start": [], "end":[]}
    for key in abq_keywords:
        for key_idx in keyword_idxs[key]:
            start_data = key_idx[0]
            idx_entries["start"].append(start_data)
    idx_entries["start"].sort()
    idx_entries["start"].append(-1)

    # update last index
    for key in abq_keywords:
        for key_idx in keyword_idxs[key]:
            id = idx_entries["start"].index(key_idx[0])
            key_idx[1] = idx_entries["start"][id+1]

    # store strings before last borderline
    abaqus_case_settings["header"] = data[:borderline_idxs[-1][0]]

    # store strings for each abaqus keyword
    for key in abq_keywords:
        for key_idx in keyword_idxs[key]:
            abaqus_case_settings[key] += coco_step_data[key_idx[0]:key_idx[1]]

    return abaqus_case_settings

def coupling_interface_template(interface_name):
    temp = [
        f"** Name: {interface_name}-PRESSURE   Type: Pressure",
        "*Dsload, op=NEW",
        f"{interface_name}, PNU, 1.",
        f"** Name: {interface_name}-SURFACE_TRACTION   Type: Surface traction",
        "*Dsload, op=NEW",
        f"{interface_name}, TRVECNU, 1., 1., 0., 0.\n",
    ]
    return "\n".join(temp)

def define_coupling_loads(abq_dict, inp_settings):
    # interface name
    interface_name = inp_settings["mp_in"][0]

    # flag is true if loads list is empty
    are_loads_empty = len(abq_dict["loads"]) == 0

    # flag is true if interface name is found in loads description
    are_coupling_loads_defined = re.search(r"[\s]+" + interface_name, abq_dict["loads"])

    load_header = []
    if are_loads_empty:
        load_header = "** LOADS\n**\n"

    load_content = []
    if not are_coupling_loads_defined:
        load_content = coupling_interface_template(interface_name)
        load_content += "**\n"

    coupled_load_setting = load_header + load_content

    abq_dict["loads"] += "".join(coupled_load_setting)

def modify_step_settings(abq_dict, inp_settings):
    check_line = abq_dict["step"]

    newline = '\n'
    lines = check_line.replace(newline, ';').split(';')
    t_new = str(inp_settings["delta_time"])
    inc_new = f"inc={str(inp_settings["inc"])}"

    # find *Step line
    step_re = r"\*(Step).*"
    for id, line in enumerate(lines):
        if re.search(step_re, line):
            step_id = id
            break
    times_id = step_id+2

    # modify step inp_settings
    step_vals = [line for line in lines[step_id].split(',')]
    step_vals[1] = ' name=coconut_abaqusStep'
    inc_re = r'inc'
    if re.search(inc_re, step_vals[-1]):
        step_vals[-1] = inc_new
    else:
        step_vals.append(inc_new)
    lines[step_id] = ",".join(step_vals)

    # modify time inp_settings
    time_setting_vals = [line for line in lines[times_id].split(',')]
    time_setting_vals[0] = t_new
    time_setting_vals[1] = t_new
    lines[times_id] = ",".join(time_setting_vals)

    # join lines
    new_lines = newline.join(lines)
    abq_dict["step"] = new_lines

def iter_write_abq_dict(f, abq_dict, print_interaction=True):
    # for step in abq_dict["step"]:
    f.write(abq_dict["step"])
    # for bcs in abq_dict["boundary conditions"]:
    f.write(abq_dict["boundary conditions"])

    # should interaction be printed
    if print_interaction:
        # for interaction in abq_dict["interactions"]:
        f.write(abq_dict["interactions"])

    # for load in abq_dict["loads"]:
    f.write(abq_dict["loads"])
    # for control in abq_dict["controls"]:
    f.write(abq_dict["controls"])
    # for output in abq_dict["output requests"]:
    f.write(abq_dict["output requests"])
    # for field in abq_dict["field output"]:
    f.write(abq_dict["field output"])
    # for history in abq_dict["history output"]:
    f.write(abq_dict["history output"])

def write_start_inp(file_out, abq_dict):
    with open(file_out, "w") as f:
        f.write(abq_dict["header"])
        f.write(f"** ----------------------------------------------------------------\n")
        f.write(f"**\n")
        iter_write_abq_dict(f, abq_dict)

def write_restart_inp(file_out, abq_dict):
    with open(file_out, "w") as f:
        f.write(f"*HEADING\n")
        f.write(f"*RESTART, READ\n")
        f.write(f"** ----------------------------------------------------------------\n")
        f.write(f"**\n")
        iter_write_abq_dict(f, abq_dict, print_interaction=False)

def preprocess_inp_file(inp_settings, case_inp_path, start_inp_path, restart_inp_path):
    abq_dict = read_input_file(case_inp_path)

    # insert coupling loads
    define_coupling_loads(abq_dict, inp_settings)

    # modify
    modify_step_settings(abq_dict, inp_settings)

    # write start file
    write_start_inp(start_inp_path, abq_dict)

    # modify
    modify_step_settings(abq_dict, inp_settings)

    # write restart file
    write_restart_inp(restart_inp_path, abq_dict)