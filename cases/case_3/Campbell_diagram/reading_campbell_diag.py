import os, sys
import numpy as np
import matplotlib.pyplot as plt
import itertools
import pickle

# making the path relative to the project
local_folder = os.getcwd()[:]
os.chdir('..\..\..')
project_folder = os.getcwd()[:]
os.chdir(local_folder)
sys.path.append(project_folder)

from wrappers.ansys_apdl_wrapper import pyAPDL



def read_freq(rot_speed, prep_file):
    ''' Function to send rotational speed and get the frequencies
    '''
    var_dict['rotational_speed'] = rot_speed # rotational speed
    folder = os.path.join(local_folder,'simulations_1') # folder to save the simulations

    apdl.update_input_variables(var_dict)
    f = apdl.read_frequencies(result_folder=result_folder)
    return f



# setting some important path
local_folder = os.getcwd()   # getting local path
filename = 'base_apdl_script.dat' # Ansys APDL base script name
filepath = os.path.join(local_folder,filename)  # set the location of ansys APDL base script
temp_folder_prefix = r'simulations_2\rot_' # prefix folder temporary folders with diferrent rotation speeds
prep_file = 'real_sector_course_prep'

apdl = pyAPDL(filepath,temp_folder_prefix)
apdl.set_ansys_apdl_path(r'C:\Program Files\ANSYS Inc\v181\ansys\bin\winx64\ansys181.exe') # location of Ansys.exe


var_dict = {}
var_dict['apdl_pre_file'] = prep_file   # APDL script with mesh and boundary conditions
var_dict['file_directory'] = local_folder 
var_dict['number_of_harmonics'] = 40 # number of max harmonic index, starts from 0
var_dict['number_of_modes'] = 12 # number of modes to compute

# evaluating different speeds
#rot_list = [0.1].extend(list(np.arange(100,5001,100)))
#rot_list = np.arange(100,5000,100)
rot_list = [0.1]
rot_list.extend(list(np.arange(100,2001,50)))
freq_list = []
freq_dict = {}

for rot in rot_list:
    result_folder = os.path.join(local_folder ,temp_folder_prefix + str(rot))
    freq = read_freq(rot, prep_file)
    for harm_id, freq_list in apdl.freq_dict.items():
        if  harm_id in  freq_dict:
            freq_dict[harm_id].append(freq_list)
        else:
            freq_dict[harm_id] = [freq_list]

ax = plt.axes()

def create_color_list():
    color_list = []
    for color in itertools.product([0,0.25,0.5,0.75,1], repeat=3):
        color_list.append(color)
    return color_list 

color_list = create_color_list()


def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

# sample usage
save_object(freq_dict, 'freq_dict_2.pkl')

for harm_id, freq_list in freq_dict.items():
    bool_label = True
    for freq in np.array(freq_list).T:
        if bool_label:
            ax.plot(rot_list[0:10],freq[0:10],'-o', color = color_list[harm_id], label='Nodal Diam. ' + str(harm_id))
            bool_label = False
        else:    
             ax.plot(rot_list[0:10],freq[0:10],'-o', color = color_list[harm_id])
    if harm_id>5:
        break

ax.legend()
plt.ylabel('Natural Frequency')        
plt.xlabel('Rotation speed')
plt.title('Campbell diagram')
plt.show()
