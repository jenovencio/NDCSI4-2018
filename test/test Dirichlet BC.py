import amfe
import matplotlib.pyplot as plt
import pandas as pd
import scipy.linalg as linalg
from scipy.sparse import linalg as splinalg
import numpy as np

from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
#%matplotlib notebook

mshfile = amfe.amfe_dir('meshes/test_meshes/2_partitions_quad_mesh.msh')
#mshfile = amfe.amfe_dir('meshes/test_meshes/3_partition_2d_blade_quad_mesh.msh')

m = amfe.Mesh()
m.import_msh(mshfile)


my_comp = amfe.CraigBamptonComponent()
my_comp.set_mesh_obj(m)

my_material = amfe.KirchhoffMaterial(E=210E9, nu=0.3, rho=7.86E3, plane_stress=True, thickness=1.0)

my_comp.set_domain(3,my_material)

K, f = my_comp.assembly_class.assemble_k_and_f()
K_, f = my_comp.assembly_class.assemble_k_and_f_neumann()
M = my_comp.assembly_class.assemble_m()




def get_dirichlet_dofs(submesh_obj,direction ='xyz',id_matrix=None):
    
    x_dir = 0
    y_dir = 1
    z_dir = 2
    
    dofs_to_keep = []
    if 'x' in direction:
        dofs_to_keep.append(x_dir)

    if 'y' in direction:
        dofs_to_keep.append(y_dir)
    
    if 'z' in direction:
        dofs_to_keep.append(z_dir)
    
    
    
    dir_nodes = submesh_obj.global_node_list
    
    dir_dofs = []
    for node, dofs in id_matrix.items():
        if node in dir_nodes:
            local_dofs = []
            for i in dofs_to_keep:
                try:
                    local_dofs.append(dofs[i])
                except:
                    print('It is not possible to issert dof %i as dirichlet dof' %i)
            dir_dofs.extend(local_dofs)
    return dir_dofs

dirsub = m.get_submesh('phys_group', 1)
dirsub2 = m.get_submesh('phys_group', 2)

id_matrix = my_comp.assembly_class.id_matrix
dir_dofs = amfe.get_dirichlet_dofs(dirsub, direction ='xyz', id_matrix=id_matrix)
dir_dofs_2 = amfe.get_dirichlet_dofs(dirsub2, direction ='x', id_matrix=id_matrix)

K, M, f = my_comp.insert_dirichlet_boundary_cond(K,M,f,dir_dofs, value = 0.0)


K, M, f = my_comp.insert_dirichlet_boundary_cond(K,M,f,dir_dofs_2, value = 1.0)

print(my_comp.check_symmetry(K))

u = splinalg.spsolve(K,f)



print(u[dir_dofs])
print(u[dir_dofs_2])

teste = 1