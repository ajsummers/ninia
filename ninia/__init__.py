# Copyright 2022 Alex Summers
# See LICENSE file for more information

import sys

if sys.version_info[0] == 2:
    raise ImportError('Ninia requires Python 3.10. This is Python 2.')

__all__ = ['Relax', 'Control', 'System', 'Electrons', 'Cell', 'Ions',
           'Job', 'parse_sisso_eqn', 'gen_sisso', 'run_sisso', 'eV_Ry',
           'Bohr_A']

__version__ = '0.1.4'
__author__ = 'Alex Summers'

from ninia.relax import Relax
from ninia.sisso import gen_sisso, run_sisso
from ninia.utils import parse_sisso_eqn, Control, System, Electrons, Cell, Ions, Job
from ninia.constants import eV_Ry, Bohr_A

