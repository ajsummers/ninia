
from dataclasses import dataclass
from jinja2 import Environment, FileSystemLoader

environment = Environment(loader=FileSystemLoader('ninia/input/'))
template = environment.get_template('relax.jinja2')






content = template.render(
    control=control,
    system=system,
    electrons=electrons,
    atomic_species=atomic_species,
    cell_parameters=cell_parameters,
    atomic_positions=atomic_positions,
    k_points=k_points,
)
