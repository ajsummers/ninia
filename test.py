from ase.build import molecule, add_adsorbate, hcp0001
from ase.visualize import view
import numpy as np

surface = hcp0001('Ru', size=(4, 2, 4), a=2.7059, c=4.2815)
ad = molecule('NH2')
ad.rotate(180, 'x')
add_adsorbate(surface, ad, 2.0, 'hcp')


def get_atomic_info(atom_obj):

    atomic_positions = ''
    positions = atom_obj.get_positions().tolist()
    symbols = atom_obj.get_chemical_symbols()
    unique_symbols = list(set(symbols))
    atom_count = len(positions)
    for atom_set in zip(symbols, positions):
        atomic_positions += f'   {atom_set[0]}\t{np.round(atom_set[1][0], 8):.8f}\t{np.round(atom_set[1][1], 8):.8f}'
        atomic_positions += f'\t{np.round(atom_set[1][2], 8):.8f}\n'
    return atomic_positions, unique_symbols, atom_count


atomic_info = get_atomic_info(surface)
print(atomic_info[0])
print(atomic_info[1:])

prefix = 'QE_ninia_test'
output_dir = '/home/ajs0201/workQE/output/'
pseudo_dir = '/home/ajs0201/workQE/pseudo/'

num_atoms = atomic_info[2]
num_elem = len(atomic_info[1])
ecutwfc = 30.0
input_dft = 'beef'

conv_thr = '1.0d-8'
mixing_beta = '0.15d0'


