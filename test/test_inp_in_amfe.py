# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 17:58:15 2018

@author: ge72tih
"""

# importing lib and setting a list of meshes to be tested
import sys 
sys.path.append(r'C:\AMfe')
sys.path.append(r'C:\Projects')


import amfe
from wrappers.read_abaqus_mesh import *

filepath = r'C:\Projects\data\beampl.inp'
nodes_dict, elem_list, nset_list  = read_inp(filepath)
el_df, node_idx = create_amfe_elem_data_frame(elem_list)
nodes = create_amfe_node_array(nodes_dict)

m2 = amfe.Mesh()
m2.el_df = el_df
m2.node_idx = node_idx
m2.nodes = nodes
m2.no_of_dofs_per_node = 3
m2._update_mesh_props()

# creating a mechanical component
my_comp = amfe.MechanicalSystem()

# creating material
my_material = amfe.KirchhoffMaterial(E=210E9, nu=0.3, rho=7.86E3,  plane_stress=False)

# setting mesh object and selecting problem domain
my_comp.set_mesh_obj(m2)
domain = my_comp.set_domain('B1',my_material)

my_comp.assembly_class.compute_element_indices()
K, f = my_comp.assembly_class.assemble_k_and_f()