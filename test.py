
from jinja2 import Environment, BaseLoader
import pkg_resources
from ninia.utils import Control, Job

input_string = pkg_resources.resource_string(__name__, 'ninia/input/slurm.jinja2')
input_template = Environment(loader=BaseLoader).from_string(input_string.decode('utf-8'))
string = input_template.render(job=Job(), control=Control())
print(string)
