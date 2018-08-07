import os, sys
import numpy as np

# making the path relative to the project
local_folder = os.getcwd()[:]
os.chdir('..\..\..')
project_folder = os.getcwd()[:]
os.chdir(local_folder)
sys.path.append(project_folder)

from wrappers.ansys_apdl_wrapper import pyAPDL



def calc_freq(rot_speed, prep_file):
    ''' Function to send rotational speed and get the frequencies
    '''
    var_dict['rotational_speed'] = rot_speed # rotational speed

    folder = os.path.join(local_folder,'simulations_2') # folder to save the simulations
    apdl.update_input_variables(var_dict)
    apdl.write_apdl_script(folder,rot_speed)

    
    apdl.run_simulation()
    return None




# seeting some important path
local_folder = os.getcwd()   # getting local path
filename = 'base_apdl_script.dat' # Ansys APDL base script name
filepath = os.path.join(local_folder,filename)  # set the location of ansys APDL base script
temp_folder_prefix = 'rot_' # prefix folder temporary folders with diferrent rotation speeds

prep_file = 'real_sector_course_prep'


apdl = pyAPDL(filepath,temp_folder_prefix)
apdl.set_ansys_apdl_path(r'C:\Program Files\ANSYS Inc\v181\ansys\bin\winx64\ansys181.exe') # location of Ansys.exe


var_dict = {}
var_dict['apdl_pre_file'] = prep_file   # APDL script with mesh and boundary conditions
var_dict['file_directory'] = local_folder 
var_dict['number_of_harmonics'] = 40 # number of max harmonic index, starts from 0
var_dict['number_of_modes'] = 12 # number of modes to compute

# evaluating different speeds
#rot_list = [0.1]
#rot_list.extend(list(np.arange(100,2001,50)))
rot_list= list(np.arange(1500,2001,50))

for rot in rot_list:
    freq = calc_freq(rot, prep_file)



