# Ninia

A small Python wrapper for setting up Quantum Espresso input files. More functionality may be added later.

<p>Currently, there is only assumed support for hexagonal close packed (HCP) crystal structures. 
Support for other structures may be tested but *should not* be expected.</p>

---
#### Example usage:

```python
# Import necessary modules:
from ase.build import molecule, add_adsorbate, hcp0001
from ase.visualize import view

# Set up geometry using ASE:
surface = hcp0001('Ru', size=(4, 2, 4), a=2.7059, c=4.2815)
ad = molecule('NH2')
ad.rotate(180, 'x')
add_adsorbate(surface, ad, 2.0, 'hcp')

view(surface, viewer='x3d')  # Specific viewer for use in Jupyter
```
<p>This will display a view of the geometry we have created. More information
about ASE (Atomic Simulation Environment) can be found at their homepage:
<a href="https://wiki.fysik.dtu.dk/ase/">https://wiki.fysik.dtu.dk/ase/</a></p>

Then you can start using ninia to convert this geometry into an input file:
```python
from ninia import relax
calc = relax.Relax(prefix='Ru_test', functional='beef',
                   pseudodir='/home/ajs0201/workQE/pseudo')
calc.load_geometry(surface)
calc.set_directories(outputdir='/home/ajs0201/workQE/output')
# Ninia assumes the current script directory as the input directory
# if none is given.
calc.set_parameters(mixing_beta=0.15)
```
