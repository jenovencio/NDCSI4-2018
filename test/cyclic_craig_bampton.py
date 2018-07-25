# importing lib and setting a list of meshes to be tested
import sys 
import amfe
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import linalg as splinalg
import copy


# making the path relative to the project
local_folder = os.getcwd()[:]
os.chdir('..')
project_folder = os.getcwd()[:]
os.chdir(local_folder)
sys.path.append(r'project_folder')

mesh_inp = os.path.join(project_folder,r'data\ansys_simple_blade_disc_files\simple_blade_disc.inp')


m = amfe.Mesh()
m.import_inp(mesh_inp,1.0)



my_comp = amfe.CraigBamptonComponent()
my_comp.set_mesh_obj(m)

my_material = amfe.KirchhoffMaterial(E=210E9, nu=0.3, rho=7.86E3, plane_stress=True, thickness=1.0)
my_comp.set_domain('SOLID_1_1_SOLID_ELSET',my_material)


#-----------------------------------------------------------------------------------------------------
if False:
    K, f = my_comp.assembly_class.assemble_k_and_f()
    K_, f = my_comp.assembly_class.assemble_k_and_f_neumann()
    M = my_comp.assembly_class.assemble_m()
    amfe.save_object(K,'K.pkl')
    amfe.save_object(M,'M.pkl')
    amfe.save_object(f, 'f.pkl')
    amfe.save_object(my_comp,'my_comp.pkl')
else:
    K = amfe.load_obj('K.pkl')
    M = amfe.load_obj('M.pkl')
    f = amfe.load_obj('f.pkl')
    my_comp = amfe.load_obj('my_comp.pkl')


dirsub = m.get_submesh('phys_group', 'DIRICHLET_ELSET')
cyclic_low = m.get_submesh('phys_group', 'LOW_ELSET')
cyclic_high = m.get_submesh('phys_group', 'HIGH_ELSET')
cyclic_top_low = m.get_submesh('phys_group', 'TOP_LOW_ELSET')
cyclic_top_high = m.get_submesh('phys_group', 'TOP_HIGH_ELSET')

id_matrix = my_comp.assembly_class.id_matrix
dir_dofs = amfe.get_dirichlet_dofs(dirsub, direction ='xyz', id_matrix=id_matrix)
low_dofs = amfe.get_dirichlet_dofs(cyclic_low, direction ='xyz', id_matrix=id_matrix)
high_dofs = amfe.get_dirichlet_dofs(cyclic_high, direction ='xyz', id_matrix=id_matrix)
top_low_dofs = amfe.get_dirichlet_dofs(cyclic_top_low, direction ='xyz', id_matrix=id_matrix)
top_high_dofs = amfe.get_dirichlet_dofs(cyclic_top_high, direction ='xyz', id_matrix=id_matrix)


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

#--------------------------------------------------------------------------
# compute the cyclic Craig Bampton reduction
K, M, f = my_comp.insert_dirichlet_boundary_cond(K,M,f,dir_dofs)


Tcg = my_comp.compute_cyclic_Craig_Bampton_Red(M, K, master_dofs, low_dofs, high_dofs, num_of_sector, no_of_modes=6, harm_index=0)


M_red = Tcg.T.dot(M).dot(Tcg)
K_red = Tcg.T.dot(K).dot(Tcg)
C_red = Tcg.T.dot(0.003*K).dot(Tcg)

def writeMatrix(M,filename):

    ndof, garbage = M.shape
    with open(filename,'w') as f:

        f.write(str(ndof)) 
        f.write('\n') 
        for columns in M.todense().T:
            for i in columns.tolist()[0]:
                f.write(str(i))
                f.write('\n') 


M_file = 'M_red_0.dat'
K_file = 'K_red_0.dat'
C_file = 'C_red_0.dat'

writeMatrix(M_red,M_file)
writeMatrix(K_red,K_file)
writeMatrix(C_red,C_file)