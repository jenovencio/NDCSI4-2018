# importing lib and setting a list of meshes to be tested
import sys 
import amfe
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import linalg as splinalg
import copy
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import mpl_toolkits.mplot3d as a3

# making the path relative to the project
local_folder = os.getcwd()[:]
os.chdir('..\..')
project_folder = os.getcwd()[:]
os.chdir(local_folder)
sys.path.append(project_folder)


from utils.write import *

mesh_inp = os.path.join(project_folder,r'cases\case_1\simple_blade_disc.inp')


m = amfe.Mesh()
m.import_inp(mesh_inp,1.0)

my_comp = amfe.CraigBamptonComponent()
my_comp.set_mesh_obj(m)

my_material = amfe.KirchhoffMaterial(E=210E9, nu=0.3, rho=7.86E3, plane_stress=True, thickness=1.0)
my_comp.set_domain('SOLID_1_1_SOLID_ELSET',my_material)


#-----------------------------------------------------------------------------------------------------
# Assembly the system operator

try:
    K = amfe.load_obj('K.pkl')
    M = amfe.load_obj('M.pkl')
    f = amfe.load_obj('f.pkl')
    my_comp = amfe.load_obj('my_comp.pkl')
except:
    K, f = my_comp.assembly_class.assemble_k_and_f()
    K_, f = my_comp.assembly_class.assemble_k_and_f_neumann()
    M = my_comp.assembly_class.assemble_m()
    amfe.save_object(K,'K.pkl')
    amfe.save_object(M,'M.pkl')
    amfe.save_object(f, 'f.pkl')
    amfe.save_object(my_comp,'my_comp.pkl')


#dirsub = m.get_submesh('phys_group', 'DIRICHLET_ELSET')
#cyclic_low = m.get_submesh('phys_group', 'LOW_ELSET')
#cyclic_high = m.get_submesh('phys_group', 'HIGH_ELSET')
cyclic_top_low = m.get_submesh('phys_group', 'TOP_LOW_ELSET')
cyclic_top_high = m.get_submesh('phys_group', 'TOP_HIGH_ELSET')


#id_matrix = my_comp.assembly_class.id_matrix
#dir_dofs = amfe.get_dirichlet_dofs(dirsub, direction ='xyz', id_matrix=id_matrix)
#low_dofs = amfe.get_dirichlet_dofs(cyclic_low, direction ='xyz', id_matrix=id_matrix)
#high_dofs = amfe.get_dirichlet_dofs(cyclic_high, direction ='xyz', id_matrix=id_matrix)
#top_low_dofs = amfe.get_dirichlet_dofs(cyclic_top_low, direction ='xyz', id_matrix=id_matrix)
#top_high_dofs = amfe.get_dirichlet_dofs(cyclic_top_high, direction ='xyz', id_matrix=id_matrix)


dir_dofs = my_comp.get_dofs('DIRICHLET_ELSET')
low_dofs = my_comp.get_dofs('LOW_ELSET')
high_dofs = my_comp.get_dofs('HIGH_ELSET')
top_low_dofs = my_comp.get_dofs('TOP_LOW_ELSET')
top_high_dofs = my_comp.get_dofs('TOP_HIGH_ELSET')


sector_angle = 45
contact= amfe.Cyclic_Contact(cyclic_top_low, cyclic_top_high, sector_angle = sector_angle)
write_contact('contact_setup.dat',contact)

# concatenating cyclic dofs
all_low_dofs =[]
all_low_dofs.extend(low_dofs)
all_low_dofs.extend(top_low_dofs)
all_high_dofs = []
all_high_dofs.extend(high_dofs)
all_high_dofs.extend(top_high_dofs)


# remove Dirichle dofs at cyclic dofs 
dir_dofs = list(set(dir_dofs).difference(all_low_dofs))
dir_dofs = list(set(dir_dofs).difference(all_high_dofs))
num_of_sector = 8

master_dofs = []
master_dofs.extend(top_low_dofs)
master_dofs.extend(top_high_dofs)
master_dofs.extend(dir_dofs)


harm_id_list = [0]
no_of_modes = 10
alpha = 0.003
beta = 0.0

for harm_id in harm_id_list:
    # Building the Cyclic Reduced System for each Harminc
    Tcg, Tcg_inv  = my_comp.compute_cyclic_Craig_Bampton_Red(M, K, master_dofs, low_dofs, high_dofs, num_of_sector, no_of_modes=no_of_modes, harm_index=harm_id)


    # building maps to reduced system and global system
    inv_id_matrix = my_comp.create_inv_id_matrix()
    local2global_master = my_comp.local2global_master 



    filename = 'master_map.dat'
    if not os.path.isfile(filename):
        write_map_matrix(filename,inv_id_matrix,local2global_master,label=False)

    # create K and M matrix with Dirichlet B.C
    K, M, f = my_comp.insert_dirichlet_boundary_cond(K,M,f,dir_dofs)
    M_red = Tcg.T.dot(M).dot(Tcg)
    K_red = Tcg.T.dot(K).dot(Tcg)
    C_red = Tcg.T.dot(alpha*K + beta*M).dot(Tcg)

    

    real_tag = '_red_' + str(harm_id) + '_real.dat'
    imag_tag = '_red_' + str(harm_id) + '_imag.dat'

    M_file = 'M' + real_tag
    K_file = 'K' + real_tag
    C_file = 'C' + real_tag

    Mi_file = 'M' + imag_tag
    Ki_file = 'K' + imag_tag
    Ci_file = 'C' + imag_tag

    writeMatrix(np.real(M_red),M_file)
    writeMatrix(np.real(K_red),K_file)
    writeMatrix(np.real(C_red),C_file)

    writeMatrix(np.imag(M_red),Mi_file)
    writeMatrix(np.imag(K_red),Ki_file)
    writeMatrix(np.imag(C_red),Ci_file)


n, m1 = Tcg.shape
with open('case_setup.dat','w') as f:
    
    num_of_lines = 3 + len(harm_id_list)
    f.write(str(num_of_lines) +'\n')
    f.write(str(m1) +'\n')
    
    f.write(str(np.deg2rad(sector_angle)) +'\n')
    for harm_id in harm_id_list:
        f.write(str(harm_id) +'\n')



