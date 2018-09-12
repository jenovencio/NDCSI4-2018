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


dirsub = m.get_submesh('phys_group', 1)
dirsub2 = m.get_submesh('phys_group', 2)

id_matrix = my_comp.assembly_class.id_matrix
dir_dofs = amfe.get_dirichlet_dofs(dirsub, direction ='xyz', id_matrix=id_matrix)
dir_dofs_2 = amfe.get_dirichlet_dofs(dirsub2, direction ='x', id_matrix=id_matrix)

value_1 = 0.0
value_2 = 2.0

K, M, f = my_comp.insert_dirichlet_boundary_cond(K,M,f,dir_dofs, value = value_1)
K, M, f = my_comp.insert_dirichlet_boundary_cond(K,M,f,dir_dofs_2, value = value_2)

Kii, S, S_inv, T, fi = my_comp.create_selection_operator(dir_dofs, K, f, remove = True)

ui = splinalg.spsolve(Kii,fi)

u = T(ui)



u_dir_1 = u[dir_dofs]
u_dir_2 = u[dir_dofs_2]
tol = 1.0e-8
for ui in u_dir_1:
    if abs(ui - value_1)>tol:
        raise('Error in the selection operator method')

for ui in u_dir_2:
    if abs(ui - value_2)>tol:
        raise('Error in the selection operator method')