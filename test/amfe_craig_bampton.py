import amfe
import matplotlib.pyplot as plt
import pandas as pd
import scipy.linalg as linalg

mshfile = amfe.amfe_dir('meshes/test_meshes/2_partitions_quad_mesh.msh')

m = amfe.Mesh()
m.import_msh(mshfile)

my_comp = amfe.CraigBamptonComponent()
my_comp.set_mesh_obj(m)

my_material = amfe.KirchhoffMaterial(E=210E9, nu=0.3, rho=7.86E3, plane_stress=True, thickness=1.0)
#domain = m.set_domain('phys_group', 3)
#domain.set_material(my_material)


my_comp.set_domain(3,my_material)

my_comp.apply_neumann_boundaries(2,1e8, 'normal')
K, f = my_comp.assembly_class.assemble_k_and_f()
K_, f = my_comp.assembly_class.assemble_k_and_f_neumann()
M = my_comp.assembly_class.assemble_m()

dirsub = m.get_submesh('phys_group', 1)
dir_nodes = dirsub.global_node_list
id_matrix = my_comp.assembly_class.id_matrix
dir_dofs = []
for node, dof in id_matrix.items():
    if node in dir_nodes:
        dir_dofs.extend(dof)

        
master_sub = m.get_submesh('phys_group', 2)
master_nodes = master_sub.global_node_list
id_matrix = my_comp.assembly_class.id_matrix
master_dofs = []
slave_dofs = []
for node, dof in id_matrix.items():
    if node in master_nodes:
        master_dofs.extend(dof)        
    else:
        slave_dofs.extend(dof)

slave_dofs = list(set(slave_dofs).difference(dir_dofs))
master_dofs.extend(dir_dofs)
    
K, M, f = my_comp.insert_dirichlet_boundary_cond(K,M,f,dir_dofs)

no_of_modes = 20
T, T_local, P, K_local, M_local = my_comp.compute(M, K, master_dofs, slave_dofs, no_of_modes=no_of_modes)


Mr = T.T.dot(M.todense()).dot(T)
Kr = T.T.dot(K.todense()).dot(T)
ff = T.T.dot(f)

no_of_modes = 10
omega, V_dynamic = linalg.eig(K.todense(), M.todense())
omega_red, V_dynamic_red = linalg.eig(Kr, Mr)

omega = np.sort(omega)
omega_red = np.sort(omega_red)

plt.plot(omega_red, 'o')
plt.plot(omega, 'x')
plt.show()