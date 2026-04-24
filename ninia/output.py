
from ase import Atoms
import numpy as np
import glob
import re


def get_final_atoms(output_file, input_file=None):

    files = glob.glob(output_file)
    if files:
        with open(files[0], 'r') as f:
            handle = f.read()
    else:
        raise IOError('No output file found')

    match = r'(?<=ATOMIC_POSITIONS)(.*?)(?=(End final coordinates|\n\n))'
    matches = re.findall(match, handle, flags=re.S)

    element_list, x_list, y_list, z_list = [[], [], [], []]
    last_match = matches[-1]

    formatted = '\n'.join(last_match[0].split('\n')[1:])
    num_atoms = len(formatted.split('\n')) - 1
    for line in formatted.split('\n')[:-1]:
        line_list = line.split()
        element_list.append(line_list[0])
        x_list.append(float(line_list[1]))
        y_list.append(float(line_list[2]))
        z_list.append(float(line_list[3]))

    position_list = [tuple(group) for group in zip(x_list, y_list, z_list)]
    atoms = Atoms(element_list, positions=position_list)

    if input_file is not None:

        with open(input_file, 'r') as f:
            handle = f.read()

        pattern = (
            r'CELL_PARAMETERS angstrom\s*'
            r'((?:[ \t]*[-+]?[\d.]+(?:[eE][-+]?\d+)?\s+'
            r'[-+]?[\d.]+(?:[eE][-+]?\d+)?\s+'
            r'[-+]?[\d.]+(?:[eE][-+]?\d+)?\n){3})'
        )

        # Find all occurrences of the cell parameters block.
        matches = re.findall(pattern, handle)
        if matches:
            # Get the final occurrence.
            final_block = matches[-1]
            # Split the block into individual lines and convert each to a list of floats.
            cell_lines = final_block.strip().splitlines()
            cell_parameters = [list(map(float, line.split())) for line in cell_lines]
        else:
            raise AttributeError('No cell parameters found')

        cell_parameters = np.array(cell_parameters)

        atoms.set_cell(cell_parameters)

    return atoms


def get_final_energy(output_file):

    files = glob.glob(output_file)
    if files:
        with open(files[0]) as f:
            content = f.read()
        match = re.search(r'Final energy\s*=\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*Ry', content)
        return float(match.group(1)) if match else np.nan
    return np.nan

