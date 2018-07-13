# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 12:45:00 2018

@author: SIV0021
"""

import sys
sys.path.append(r'D:\NDCSI\framework')


import numpy as np


""" Read the dof array  """

dof_map_str = open(r'D:\NDCSI\framework\data\calculix_files_exportK_readinpython\without_prestress\stdrectblock_quadface_hexmesh_withprestressforce.dof',"r")
dof_map_lines = dof_map_str.readlines()

dof_map = []
for line in dof_map_lines:
    value = line.split(".");
    dof_map.append([ int(value[0]),int(value[1]) ])
    
    
""" read the mass matrix """
    
mass_mat_str =   open(r'D:\NDCSI\framework\data\calculix_files_exportK_readinpython\without_prestress\stdrectblock_quadface_hexmesh_withprestressforce.mas',"r")    
mass_mat_lines = mass_mat_str.readlines()     
#for i in range(0,len(mass_mat_lines)):

#    tmp = str.split(mass_mat_lines[i])
#    mass_mat_list[i,:] = [int(tmp[0]),int(tmp[1]),float(tmp[2])]
mass_mat = []
n = len(mass_mat_lines) 
for line in mass_mat_lines:
    value = line.split()
    mass_mat.append([int(value[0]),int(value[1]),float(value[2])])
    
    
    
""" read the stiffness matrix """

stiff_mat_str =   open(r'D:\NDCSI\framework\data\calculix_files_exportK_readinpython\without_prestress\stdrectblock_quadface_hexmesh_withprestressforce.sti',"r")    
stiff_mat_lines = stiff_mat_str.readlines()     
#for i in range(0,len(mass_mat_lines)):

#    tmp = str.split(mass_mat_lines[i])
#    mass_mat_list[i,:] = [int(tmp[0]),int(tmp[1]),float(tmp[2])]
stiff_mat = []
for line in stiff_mat_lines:
    value = line.split()
    stiff_mat.append([int(value[0]),int(value[1]),float(value[2])])