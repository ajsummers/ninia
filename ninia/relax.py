"""
Alex Summers
The current function of this script is to create a Relax object to handle relevant information
and create a generic input file and bash script for use with Quantum Espresso.
For more information see the GitHub repository at https://github.com/ajsummers/ninia.
"""

from ninia.utils import Control, System, Electrons, Cell, Ions, Job
from typing import Type, Union, List, Tuple
from jinja2 import Environment, BaseLoader
from ase import Atom, Atoms
from ninia import utils
from importlib.resources import files
from pathlib import Path
import os

starting_dir = os.getcwd()

input_text = files(__package__).joinpath('input/relax.jinja2').read_text('utf-8')
input_template = Environment(loader=BaseLoader()).from_string(input_text)
slurm_text = files(__package__).joinpath('input/slurm.jinja2').read_text('utf-8')
slurm_template = Environment(loader=BaseLoader()).from_string(slurm_text)


class Relax:

    def __init__(self, control: Type[Control] = Control(), system: Type[System] = System(),
                 electrons: Type[Electrons] = Electrons(), cell: Type[Cell] = Cell(), ions: Type[Ions] = Ions(),
                 job: Type[Job] = Job(), geometry: Union[Type[Atom], Type[Atoms]] = None, input_dir: str = None,
                 k_points: Tuple[int] = (1, 1, 1, 0, 0, 0), job_preamble: str = None, job_command: str = None,):

        # Initialize class parameters

        self.control = control
        self.system = system
        self.electrons = electrons
        self.cell = cell
        self.ions = ions
        self.job = job
        self.job_preamble = job_preamble
        self.job_command = job_command
        self.geometry = geometry
        self.input_dir = input_dir
        self.k_points = k_points

        self.cell_parameters = None
        self.atomic_species = None
        self.atomic_positions = None

    def set_atomic_info(self, geometry: Union[Type[Atom], Type[Atoms]] = None) -> None:

        if (self.geometry is None) and (geometry is None):
            raise RuntimeError('Need to define geometry info (ase.Atoms object).')
        elif geometry is not None:
            self.geometry = geometry

        self.set_pseudodir()

        self.system.nat, self.system.ntyp, self.atomic_positions = utils.position(self.geometry)
        self.atomic_species = utils.species(self.geometry, pseudo_dir=self.control.pseudo_dir)
        self.cell_parameters = utils.cell_parameters(self.geometry)

    def set_pseudodir(self, pseudo_dir: str = None):

        if (self.control.pseudo_dir is None) and (pseudo_dir is None):
            raise RuntimeError('Need to define pseudo_dir.')
        elif pseudo_dir is not None:
            self.control.pseudo_dir = pseudo_dir

        if not os.path.isabs(self.control.pseudo_dir):
            self.control.pseudo_dir = os.path.join(starting_dir, self.control.pseudo_dir)

        self.check_directory_(self.control.pseudo_dir)

    def set_outdir(self, outdir: str = None):
        if (self.control.outdir is None) and (outdir is None):
            raise RuntimeError('Need to define outdir.')
        elif outdir is not None:
            self.control.outdir = outdir

        if not os.path.isabs(self.control.outdir):
            self.control.outdir = os.path.join(starting_dir, self.control.outdir)

        self.check_directory_(self.control.outdir)

    def set_inputdir(self, input_dir: str = None):
        if (self.input_dir is None) and (input_dir is None):
            self.input_dir = starting_dir
        elif input_dir is not None:
            self.input_dir = input_dir

        if not os.path.isabs(self.input_dir):
            self.input_dir = os.path.join(starting_dir, self.input_dir)

        self.check_directory_(self.input_dir)

    @staticmethod
    def check_directory_(directory):

        if not os.path.isdir(directory):
            raise NotADirectoryError(f'{directory} is not a valid directory.')

    def set_prefix_(self):

        input_file = os.path.join(self.input_dir, f'{self.control.prefix}.i')
        if not (os.path.isfile(input_file)):
            save_file = os.path.join(self.control.outdir, f'{self.control.prefix}.save')
            wfc_file = os.path.join(self.control.outdir, f'{self.control.prefix}.wfc1')
            if os.path.isdir(save_file) or os.path.isfile(wfc_file):
                self.control.prefix += '_1'

            i = 1
            save_file = os.path.join(self.control.outdir, f'{self.control.prefix}.save')
            wfc_file = os.path.join(self.control.outdir, f'{self.control.prefix}.wfc1')
            while os.path.isdir(save_file) or os.path.isfile(wfc_file):
                i += 1
                self.control.prefix = f'{self.control.prefix.rstrip("0123456789")}{i}'
                save_file = os.path.join(self.control.outdir, f'{self.control.prefix}.save')
                wfc_file = os.path.join(self.control.outdir, f'{self.control.prefix}.wfc1')

            return True

        else:

            return False

    def create_input(self, control: Type[Control] = None, system: Type[System] = None,
                     electrons: Type[Electrons] = None, cell: Type[Cell] = None, ions: Type[Ions] = None,
                     geometry: Union[Type[Atom], Type[Atoms]] = None):

        if control is not None:
            self.control = control
        if system is not None:
            self.system = system
        if electrons is not None:
            self.electrons = electrons
        if cell is not None:
            self.cell = cell
        if ions is not None:
            self.ions = ions
        if geometry is not None:
            self.geometry = geometry
            self.set_atomic_info(geometry)

        if self.atomic_positions is None:
            self.set_atomic_info(self.geometry)

        self.set_outdir()
        self.set_inputdir()

        if self.set_prefix_():

            input_content = input_template.render(control=self.control, system=self.system, electrons=self.electrons,
                                                  cell=self.cell, ions=self.ions, atomic_species=self.atomic_species,
                                                  cell_parameters=self.cell_parameters,
                                                  atomic_positions=self.atomic_positions,
                                                  k_points=self.k_points)

            input_file = os.path.join(self.input_dir, f'{self.control.prefix}.i')
            with open(input_file, 'w') as f:
                f.write(input_content)

            print(f'Created input file {input_file}')

        else:

            input_file = os.path.join(self.input_dir, f'{self.control.prefix}.i')
            print(f'Muted redundant input file {input_file}')

    def create_job(self, job: Type[Job] = None, template_file: str = None, extension: str = '.sh',
                   preamble: str = None, command: str = 'mpirun'):

        if job is not None:
            self.job = job

        if preamble is not None:
            self.job_preamble = preamble

        if command is not None:
            self.job_command = command

        self.check_directory_(self.input_dir)
        if self.job.input is None:
            self.job.input = os.path.join(self.input_dir, f'{self.control.prefix}.i')
        if self.job.output is None:
            self.job.output = os.path.join(self.input_dir, f'{self.control.prefix}.out')

        job_file = os.path.join(self.input_dir, f'{self.control.prefix}{extension}')
        if not (os.path.isfile(job_file)):

            if template_file is None:
                job_content = slurm_template.render(control=self.control, job=self.job, preamble=self.job_preamble,
                                                    command=self.job_command, )
            else:
                job_template = Environment(loader=BaseLoader()).from_string(Path(template_file).read_text())
                job_content = job_template.render(control=self.control, job=self.job, preamble=self.job_preamble,
                                                  command=self.job_command, )

            with open(job_file, 'w') as f:
                f.write(job_content)

            os.chdir(starting_dir)
            print(f'Created bash file {job_file}')

        else:

            print(f'Muted redundant bash file {job_file}')

    def lock_atoms(self, lock: Union[str, Tuple[int]] = None, which: Tuple[int] = (0, 0, 0)) -> None:

        if self.atomic_positions is None:
            self.set_atomic_info(self.geometry)

        self.atomic_positions = utils.lock_atoms(lock, which, self.atomic_positions)

