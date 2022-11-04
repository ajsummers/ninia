
from jinja2 import Environment, FileSystemLoader
from ninia.utils import Control, System, Electrons, Job
import os
from fnmatch import filter as flt

environment = Environment(loader=FileSystemLoader('ninia/input/'))
input_template = environment.get_template('relax.jinja2')
job_template = environment.get_template('slurm.jinja2')

control = Control(prefix='Ag_223_OH', outdir='/home/ajs0201/outdir', pseudo_dir='/home/pseudo_dir', nstep=100)
system = System(nat=20, ntyp=3, ecutwfc=95, ecutrho=500, input_dft='beef', degauss=0.1)
electrons = Electrons()
atomic_species = 'Ag ...'
cell_parameters = 'Cell Parameters ... 33.304434 4345. 34345'
atomic_positions = 'Ag 4.045 6.500 66.00 0 0 0'
k_points = ' 4 4 4 0 0 0'
job = Job(ntasks=32, partition='nova', mail_type=['FAIL', 'REQUEUE'], mail_user='ajs0201@auburn.edu', nk=4,
          input=f'/home/ajs0201/workQE/input/{control.prefix}.i',
          output=f'/home/ajs0201/workQE/input/{control.prefix}.out')

input_content = input_template.render(
    control=control,
    system=system,
    electrons=electrons,
    atomic_species=atomic_species,
    cell_parameters=cell_parameters,
    atomic_positions=atomic_positions,
    k_points=k_points,
)

job_content = job_template.render(
    control=control,
    job=job,
)

# print(input_content)
# print(job_content)

list_j2 = flt(os.listdir('./ninia/input'), '*.[Jj]inja2')
print(list_j2)
