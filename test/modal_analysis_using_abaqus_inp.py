# importing lib and setting a list of meshes to be tested

import sys 
import amfe
import os
import matplotlib.pyplot as plt
import numpy as np
import dill as pickle

# making the path relative to the project
local_folder = os.getcwd()[:]
os.chdir('..')
project_folder = os.getcwd()[:]
os.chdir(local_folder)
sys.path.append(r'project_folder')

mesh_inp = os.path.join(project_folder,r'data\abaqus_files\simple_sector\abaqus.inp')


m = amfe.Mesh()
m.import_inp(mesh_inp,1.0)

m.split_in_groups()

g = m.groups.keys()

sub_domain = m.groups['SOLID_1_1_SOLID_ELSET']

sub_domain.create_elem_dict()


vertice_matrix = sub_domain.get_element(0)
points_coord = sub_domain.parent_mesh.nodes

if False:
    ax = amfe.plot3Dmesh(m, ax=None, boundaries=True, alpha=0.2, color='blue', scale = 1000, plot_nodes=False)
    ax.view_init(-90, 0)
    #ax.set_aspect('equal', adjustable='box')
    xlim = ax.get_xlim()
    ax.set_ylim(xlim)
    ax.set_zlim(xlim)
    plt.axis('off')



load_pickle_obj = False
if not load_pickle_obj:
    my_comp = amfe.CraigBamptonComponent()
    my_comp.set_mesh_obj(m)

    my_material = amfe.KirchhoffMaterial(E=210E9, nu=0.3, rho=7.86E3, plane_stress=True, thickness=1.0)
    my_comp.set_domain('SOLID_1_1_SOLID_ELSET',my_material)
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

id_matrix = my_comp.assembly_class.id_matrix
dir_dofs = amfe.get_dirichlet_dofs(dirsub, direction ='xyz', id_matrix=id_matrix)
low_dofs = amfe.get_dirichlet_dofs(cyclic_low, direction ='xyz', id_matrix=id_matrix)
high_dofs = amfe.get_dirichlet_dofs(cyclic_high, direction ='xyz', id_matrix=id_matrix)


# remove Dirichle dofs at cyclic dofs 
dir_dofs = list(set(dir_dofs).difference(low_dofs))
dir_dofs = list(set(dir_dofs).difference(high_dofs))

# inserting dirichlet boundary contitions
K, M, f = my_comp.insert_dirichlet_boundary_cond(K,M,f,dir_dofs)

# removing dirichlet B.C in order to use sparse eigen solver 
Kii, S, S_inv, T, fi = my_comp.create_selection_operator(dir_dofs, K, f, remove = True)
Mii, Sm, Sm_inv, Tm = my_comp.create_selection_operator(dir_dofs, M, remove = True)




