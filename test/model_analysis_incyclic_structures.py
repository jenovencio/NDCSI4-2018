import amfe
import matplotlib.pyplot as plt
import pandas as pd
import scipy.linalg as linalg
from scipy.sparse import linalg as splinalg
import numpy as np

from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
#%matplotlib notebook

mshfile = amfe.amfe_dir('meshes/test_meshes/ring_connection.msh')


# importing the mesh
m = amfe.Mesh()
m.import_msh(mshfile)


# creating a Craig-Bampton Component
my_comp = amfe.CraigBamptonComponent()
my_comp.set_mesh_obj(m)

my_material = amfe.KirchhoffMaterial(E=210E9, nu=0.3, rho=7.86E3, plane_stress=True, thickness=1.0)
my_comp.set_domain(6,my_material)


# build the global matrix for Craig-Bampton Reduction
K, f = my_comp.assembly_class.assemble_k_and_f()
K_, f = my_comp.assembly_class.assemble_k_and_f_neumann()
M = my_comp.assembly_class.assemble_m()

# defining boundary conditions
dirsub = m.get_submesh('phys_group', 1)
dirsub_low = m.get_submesh('phys_group', 3)
dirsub_high = m.get_submesh('phys_group', 4)

id_matrix = my_comp.assembly_class.id_matrix


dir_dofs = amfe.get_dirichlet_dofs(dirsub, direction ='xy', id_matrix=id_matrix)
cylic_dofs_low = amfe.get_dirichlet_dofs(dirsub_low, direction ='xy', id_matrix=id_matrix)
cylic_dofs_high = amfe.get_dirichlet_dofs(dirsub_high, direction ='xy', id_matrix=id_matrix)

# remove Dirichle dofs at cyclic dofs 
dir_dofs = list(set(dir_dofs).difference(cylic_dofs_low))
dir_dofs = list(set(dir_dofs).difference(cylic_dofs_high))

# inserting dirichlet boundary contitions
K, M, f = my_comp.insert_dirichlet_boundary_cond(K,M,f,dir_dofs)


# inserting Cyclic Boundary contitions
K_mod, M_mod, f_mod = my_comp.insert_cyclic_symm_boundary_cond(K, M, f, low_dofs = cylic_dofs_low , high_dofs = cylic_dofs_high , theta = 0.0)    


Kii, S, S_inv, T, fi = my_comp.create_selection_operator(dir_dofs, K_mod, f_mod, remove = True)
Mii, Sm, Sm_inv, Tm = my_comp.create_selection_operator(dir_dofs, M_mod, remove = True)


num_of_modes = 30
omega, V_dynamic = splinalg.eigs(Kii, k=num_of_modes, M = Mii, which='SM')

# ordering eigenvaues, but it is complex
#indexes = np.argsort(omega)
#V_dynamic = V_dynamic [:,indexes]
#omega = omega[indexes]

V_dynamic = np.real(V_dynamic)
# baking to Augmented System
aug_mode_shapes = []
for i in range(num_of_modes):
    u_global = T(V_dynamic [:,i])
    u_dir = u_global[dir_dofs]
    if abs(u_dir).max()>0.0:
        raise('Error going back to Augmented system')
    aug_mode_shapes.append(u_global)

aug_mode_shapes = np.array(aug_mode_shapes)