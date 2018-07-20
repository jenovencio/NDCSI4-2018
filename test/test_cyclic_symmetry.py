import amfe
import matplotlib.pyplot as plt
import pandas as pd
import scipy.linalg as linalg
from scipy.sparse import linalg as splinalg
import numpy as np

from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
#%matplotlib notebook

mshfile = amfe.amfe_dir('meshes/test_meshes/ring_free.msh')


m = amfe.Mesh()
m.import_msh(mshfile)



#amfe.plot_mesh(m)
#plt.show()



my_comp = amfe.CraigBamptonComponent()
my_comp.set_mesh_obj(m)

my_material = amfe.KirchhoffMaterial(E=210E9, nu=0.3, rho=7.86E3, plane_stress=True, thickness=1.0)
my_comp.set_domain(6,my_material)
my_comp.apply_neumann_boundaries(2,1e8, 'normal')

K, f = my_comp.assembly_class.assemble_k_and_f()
K_, f = my_comp.assembly_class.assemble_k_and_f_neumann()
M = my_comp.assembly_class.assemble_m()


dirsub = m.get_submesh('phys_group', 1)
cyclic_low = m.get_submesh('phys_group', 3)
cyclic_high = m.get_submesh('phys_group', 4)


id_matrix = my_comp.assembly_class.id_matrix
dir_dofs = amfe.get_dirichlet_dofs(dirsub, direction ='xy', id_matrix=id_matrix)
low_dofs = amfe.get_dirichlet_dofs(cyclic_low, direction ='xy', id_matrix=id_matrix)
high_dofs = amfe.get_dirichlet_dofs(cyclic_high, direction ='xy', id_matrix=id_matrix)


my_comp.insert_cyclic_symm_boundary_cond(K, M, f, low_dofs = low_dofs, high_dofs = high_dofs, theta = 0.0)    