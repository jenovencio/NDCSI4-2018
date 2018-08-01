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


elem_id = 357
node_id = 331
cyclic_top_high.get_element(elem_id)
coord = cyclic_top_high.get_node_coord(node_id)
elem_coord = cyclic_top_high.get_elem_coord(elem_id)
normal_vec = cyclic_top_high.get_normal_to_element(elem_id)
elem_connect_to_node_dict = cyclic_top_high.create_elem_connect_to_node()
node_normal_vec = cyclic_top_high.get_normal_to_node( node_id, method = 'average')
node_normal_vec_1 = cyclic_top_high.get_normal_to_node( node_id, method = 'first')