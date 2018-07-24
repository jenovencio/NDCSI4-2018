# importing lib and setting a list of meshes to be tested

import sys 
import amfe
import os
import matplotlib.pyplot as plt
import numpy as np



# making the path relative to the project
local_folder = os.getcwd()[:]
os.chdir('..')
project_folder = os.getcwd()[:]
os.chdir(local_folder)
sys.path.append(r'project_folder')

mesh_inp = os.path.join(project_folder,r'data\abaqus_files\simple_sector\abaqus.inp')


m = amfe.Mesh()
m.import_inp(mesh_inp,1000.0)

m.split_in_groups()

g = m.groups.keys()

sub_domain = m.groups['SOLID_1_1_SOLID_ELSET']

sub_domain.create_elem_dict()


vertice_matrix = sub_domain.get_element(0)
points_coord = sub_domain.parent_mesh.nodes

hexa_list = [vertice_matrix[0:8]]
quad_list = amfe.get_quad_faces_from_hexa(hexa_list, surface_only=True)

vertice_matrix = quad_list
#amfe.plot_3D_polygon(points_coord*100, vertice_matrix, ax=None, alpha=0.5, color='grey', plot_nodes=True)
#plt.show()

amfe.plot3Dmesh(m, ax=None, boundaries=True, alpha=0.2, color='grey', plot_nodes=False)
plt.show()