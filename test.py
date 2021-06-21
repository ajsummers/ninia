from ase.build import molecule, add_adsorbate, hcp0001
from ase.visualize import view
import numpy as np

surface = hcp0001('Ru', size=(4, 2, 4), a=2.7059, c=4.2815)
ad = molecule('NH2')
ad.rotate(180, 'x')
add_adsorbate(surface, ad, 2.0, 'hcp')

view(surface)
