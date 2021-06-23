import pkg_resources
import warnings
import os

from fnmatch import filter as flt
import numpy as np
import re

starting_dir = os.getcwd()

mm_of_elements = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.0107, 'N': 14.0067,
              'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
              'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
              'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045,
              'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64,
              'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585,
              'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
              'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6,
              'I': 126.90447, 'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116,
              'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
              'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04,
              'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
              'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804,
              'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0254, 'Ac': 227.0278,
              'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 'Am': 243.0614,
              'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0951,
              'No': 259.1009, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271, 'Bh': 270, 'Hs': 269, 'Mt': 278,
              'Ds': 281, 'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294}


class Relax:

    def __init__(self, prefix=None, functional=None, pseudodir=None):

        if prefix is not None:
            self.prefix = prefix
        else:
            warnings.warn('Prefix not defined. Setting to "untitled".', UserWarning)

        if functional is not None:
            self.functional = functional
        else:
            warnings.warn('Functional not defined. Using BEEF-vdW by default.', UserWarning)
            self.functional = 'beef'

        if pseudodir is not None:
            self.check_directory_(pseudodir)
            self.pseudo_dir = pseudodir
        else:
            warnings.warn('Pseudopotential directory not specified. This will need to be defined '
                          'before loading geometry.', UserWarning)

        self.geometry = None
        self.output_dir = None
        self.pseudo_dir = None
        self.input_dir = None

        self.num_atoms = None
        self.num_elem = None
        self.ecutwfc = None

        self.conv_thr = None
        self.mixing_beta = None

        self.atomic_species = None
        self.cell_parameters = None
        self.atomic_positions = None
        self.k_points = None

        self.memory = None
        self.cpus = None
        self.hours = None

    def get_position_info_(self, ase_object):

        atomic_positions = ''
        positions = ase_object.get_positions().tolist()
        symbols = ase_object.get_chemical_symbols()
        unique_symbols = list(set(symbols))
        atom_count = len(positions)
        for atom_set in zip(symbols, positions):
            atomic_positions += f'   {atom_set[0]}\t{np.round(atom_set[1][0], 8):.8f}'
            atomic_positions += f'\t{np.round(atom_set[1][1], 8):.8f}\t{np.round(atom_set[1][2], 8):.8f}\n'

        self.num_atoms = atom_count
        self.num_elem = len(unique_symbols)
        self.atomic_positions = atomic_positions

        return unique_symbols

    def get_species_info_(self, species_list):

        # For some stupid reason a commented out exception was still killing the program
        # so I deleted it here.

        os.chdir(self.pseudo_dir)
        list_upf = flt(os.listdir(''), '*.[Uu][Pp][Ff]')
        species_string = ''

        for species in species_list:

            r = re.compile(rf'{species}[_|.]\S+\Z', flags=re.IGNORECASE)
            match = list(filter(r.match, list_upf))[0]
            mw_species = mm_of_elements[species]

            species_string += f'   {species}\t{mw_species}\t{match}\n'

        os.chdir(starting_dir)
        self.atomic_species = species_string

    def get_cell_parameters_(self, ase_object):

        supercell = ase_object.get_cell()
        supercell[2][2] = 2 * np.max(ase_object.get_positions().T[2])
        cell_parameters = ''

        for dimension in supercell:
            cell_parameters += f'{dimension[0]:.14f}\t{dimension[1]:.14f}\t{dimension[2]:.14f}\n'

        self.cell_parameters = cell_parameters

    @staticmethod
    def check_directory_(directory):

        if os.path.isdir(directory):
            pass
        else:
            raise NotADirectoryError(f'{directory} is not a valid directory.')

    def set_prefix_(self):

        os.chdir(self.output_dir)
        if (os.path.isdir(f'{self.prefix}.save')) or (os.path.isfile(f'{self.prefix}.wfc1')):
            self.prefix += '_1'
        else:
            pass

        i = 1
        while (os.path.isdir(f'{self.prefix}.save')) or (os.path.isfile(f'{self.prefix}.wfc1')):
            i += 1
            self.prefix = f'{self.prefix.rstrip("0123456789")}{i}'

        os.chdir(starting_dir)

    def set_directories(self, inputdir=None, outputdir=None, pseudodir=None):

        if (self.pseudo_dir is None) and (pseudodir is None):
            warnings.warn('Pseudopotential directory still not specified. Loading '
                          'geometry will likely result in an error.', UserWarning)
        elif (self.pseudo_dir is None) and (pseudodir is not None):
            self.check_directory_(pseudodir)
            self.pseudo_dir = pseudodir
        elif (self.pseudo_dir is not None) and (pseudodir is not None):
            print(f'Changing pseudodirectory to {pseudodir}.')
            self.check_directory_(pseudodir)
            self.pseudo_dir = pseudodir

        if inputdir is not None:
            self.check_directory_(inputdir)
            self.input_dir = inputdir
        else:
            self.input_dir = starting_dir

        if outputdir is not None:
            self.check_directory_(outputdir)
            self.output_dir = outputdir

    def load_geometry(self, ase_object):

        self.geometry = ase_object
        species_list = self.get_position_info_(self.geometry)
        self.get_species_info_(species_list)
        self.get_cell_parameters_(self.geometry)

    def set_parameters(self, ecutwfc=None, conv_thr=None, mixing_beta=None, k_points=None, functional=None):

        if (self.functional is None) and (functional is None):
            warnings.warn('Functional is still not specified. Creating input will likely result in an error.',
                          UserWarning)
        elif (self.functional is None) and (functional is not None):
            self.functional = functional
        elif (self.functional is not None) and (functional is not None):
            print(f'Changing functional to {functional}')

        if (self.ecutwfc is None) and (ecutwfc is None):
            self.ecutwfc = 30.0
        elif ecutwfc is not None:
            self.ecutwfc = ecutwfc

        if (self.conv_thr is None) and (conv_thr is None):
            self.conv_thr = '1.0d-8'
        elif conv_thr is not None:
            self.conv_thr = str(conv_thr).replace('e', 'd')

        if (self.mixing_beta is None) and (mixing_beta is None):
            self.mixing_beta = '0.7d0'
        elif mixing_beta is not None:
            self.mixing_beta = str(mixing_beta).replace('e', 'd')

        if (self.k_points is None) and (k_points is None):
            self.k_points = ' 4 4 4 0 0 0'
        elif k_points is not None:
            k_string = ''
            for point in k_points:
                k_string += f' {point}'
            self.k_points = k_string

    def create_input(self):

        self.set_prefix_()

        runtime_error = ''
        runtime_error += ' self.geometry' if not self.geometry else ''
        runtime_error += ' self.prefix' if not self.prefix else ''
        runtime_error += ' self.output_dir' if not self.output_dir else ''
        runtime_error += ' self.pseudo_dir' if not self.pseudo_dir else ''
        runtime_error += ' self.input_dir' if not self.input_dir else ''
        runtime_error += ' self.num_atoms' if not self.num_atoms else ''
        runtime_error += ' self.num_elem' if not self.num_elem else ''
        runtime_error += ' self.ecutwfc' if not self.ecutwfc else ''
        runtime_error += ' self.functional' if not self.functional else ''
        runtime_error += ' self.conv_thr' if not self.conv_thr else ''
        runtime_error += ' self.mixing_beta' if not self.mixing_beta else ''
        runtime_error += ' self.atomic_species' if not self.atomic_species else ''
        runtime_error += ' self.cell_parameters' if not self.cell_parameters else ''
        runtime_error += ' self.atomic_positions' if not self.atomic_positions else ''
        runtime_error += ' self.k_points' if not self.k_points else ''

        if runtime_error is not '':
            raise RuntimeError(f'Missing{runtime_error}')

        stream = pkg_resources.resource_stream(__name__, 'input/relax.i')
        with open(stream) as f:
            relax_fstring = f.read()

        compiled_fstring = compile(relax_fstring, '<fstring_from_file', 'eval')
        formatted_relax = eval(compiled_fstring)

        os.chdir(self.input_dir)
        with open(f'{self.prefix}.i', 'a') as f:
            f.write(formatted_relax)

        os.chdir(starting_dir)
        print(f'Created input file {self.prefix}.i')

    def create_bash(self, memory=None, cpus=None, hours=None):

        if memory is None:
            self.memory = 50
        else:
            self.memory = memory

        if cpus is None:
            self.cpus = 8
        else:
            self.cpus = cpus

        if hours is None:
            self.hours = 30
        else:
            self.hours = hours

        stream = pkg_resources.resource_stream(__name__, 'input/relax.sh')
        with open(stream) as f:
            relax_bash_fstring = f.read()

        compiled_bash_fstring = compile(relax_bash_fstring, '<fstring_from_file', 'eval')
        formatted_bash_relax = eval(compiled_bash_fstring)

        os.chdir(self.input_dir)
        with open(f'{self.prefix}.sh', 'a') as f:
            f.write(formatted_bash_relax)

        os.chdir(starting_dir)
        print(f'Created bash file {self.prefix}.sh')

