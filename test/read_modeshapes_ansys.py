import os, sys
import numpy as np
import copy

# seeting some important path
local_folder= copy.copy(os.getcwd()) 
os.chdir('..')
project_folder = copy.copy(os.getcwd())   # getting project path
os.chdir(local_folder)
sys.path.append(project_folder)

from wrappers.ansys_apdl_wrapper import pyAPDL

result_folder = os.path.join(project_folder, r'data\ansys_simple_blade_disc_files')

my_ansys = pyAPDL(result_folder = result_folder)


num_harmonics = 4
num_modes = 6
my_ansys.var_dict['number_of_harmonics']  = num_harmonics 
my_ansys .var_dict['number_of_modes']  = num_modes  
my_ansys.read_frequencies(filename_prefix='Frequency_harm_',ext='.txt')

freq_dict = my_ansys.freq_dict

modes_dict = my_ansys.read_mode_shapes()