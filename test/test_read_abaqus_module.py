# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 17:58:15 2018

@author: ge72tih
"""

# importing lib and setting a list of meshes to be tested
import sys 
sys.path.append(r'C:\Projects')

from wrappers import read_abaqus_mesh as w 

meshfile = r'C:\Projects\data\calculix_files_exportK_readinpython\MESH.inp'
meshfile = r'C:\Projects\data\calculix_files\stdrectblock_quadface_withhextet.inp'

nodes_dict, elem_list, nset_list, elset_list = w.read_inp(meshfile)