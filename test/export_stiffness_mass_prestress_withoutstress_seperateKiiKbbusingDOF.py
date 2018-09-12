# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 18:30:51 2018

@author: SIV0021
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 12:45:00 2018

@author: SIV0021
"""

import sys
sys.path.append(r'D:\NDCSI\framework')


import numpy as np

"""
without prestress
"""
""" Read the dof array  """

dof_map_str = open(r'D:\NDCSI\framework\data\calculix_files_exportK_readinpython\without_prestress\stdrectblock_quadface_hexmesh_withoutprestressforce.dof',"r")
dof_map_lines = dof_map_str.readlines()

g_dof_map = []
for line in dof_map_lines:
    value = line.split(".");
    g_dof_map.append([ int(value[0]),int(value[1]) ])
    
    
""" read the mass matrix """
    
mass_mat_str =   open(r'D:\NDCSI\framework\data\calculix_files_exportK_readinpython\without_prestress\stdrectblock_quadface_hexmesh_withoutprestressforce.mas',"r")    
mass_mat_lines = mass_mat_str.readlines()     
#for i in range(0,len(mass_mat_lines)):

#    tmp = str.split(mass_mat_lines[i])
#    mass_mat_list[i,:] = [int(tmp[0]),int(tmp[1]),float(tmp[2])]
g_mass_mat = []
n = len(mass_mat_lines) 
for line in mass_mat_lines:
    value = line.split()
    g_mass_mat.append([int(value[0]),int(value[1]),float(value[2])])
    
    
    
""" read the stiffness matrix """

stiff_mat_str =   open(r'D:\NDCSI\framework\data\calculix_files_exportK_readinpython\without_prestress\stdrectblock_quadface_hexmesh_withoutprestressforce.sti',"r")    
stiff_mat_lines = stiff_mat_str.readlines()     
#for i in range(0,len(mass_mat_lines)):

#    tmp = str.split(mass_mat_lines[i])
#    mass_mat_list[i,:] = [int(tmp[0]),int(tmp[1]),float(tmp[2])]
g_stiff_mat = []
for line in stiff_mat_lines:
    value = line.split()
    g_stiff_mat.append([int(value[0]),int(value[1]),float(value[2])])



"""
with prestress
"""
""" Read the dof array  """

dof_map_str = open(r'D:\NDCSI\framework\data\calculix_files_exportK_readinpython\with_prestress\stdrectblock_quadface_hexmesh_withprestressforce.dof',"r")
dof_map_lines = dof_map_str.readlines()

dof_map = []
for line in dof_map_lines:
    value = line.split(".");
    dof_map.append([ int(value[0]),int(value[1]) ])
    
    
""" read the mass matrix """
    
mass_mat_str =   open(r'D:\NDCSI\framework\data\calculix_files_exportK_readinpython\with_prestress\stdrectblock_quadface_hexmesh_withprestressforce.mas',"r")    
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

stiff_mat_str =   open(r'D:\NDCSI\framework\data\calculix_files_exportK_readinpython\with_prestress\stdrectblock_quadface_hexmesh_withprestressforce.sti',"r")    
stiff_mat_lines = stiff_mat_str.readlines()     
#for i in range(0,len(mass_mat_lines)):

#    tmp = str.split(mass_mat_lines[i])
#    mass_mat_list[i,:] = [int(tmp[0]),int(tmp[1]),float(tmp[2])]
stiff_mat = []
for line in stiff_mat_lines:
    value = line.split()
    stiff_mat.append([int(value[0]),int(value[1]),float(value[2])])
    
    
############################################################################################################################
############################################################################################################################


master_nodes = [6040,      6041,      6042,      6043,      6044,      6045,      6046,      6047,      6048,      6049,      6050,      6051,      6052,      6053,      6054,      6055,
      6056,      6057,      6058,      6059,      6060,      6061,      6062,      6063,      6064,      6065,      6066,      6067,      6068,      6069,      6070,      6071,
      6072,      6073,      6074,      6075,      6076,      6077,      6078,      6079,      6080,      6081,      6082,      6083,      6084,      6085,      6086,      6087,
      6088,      6089,      6090,      6091,      6092,      6093,      6094,      6095,      6096,      6097,      6098,      6099,      6100,      6101,      6102,      6103,
      6104,      6105,      6106,      6107,      6108,      6109,      6110,      6111,      6112,      6113,      6114,      6115,      6116,      6117,      6118,      6119,
      6120,      6121,      6122,      6123,      6124,      6125,      6126,      6127,      6128,      6129,      6130,      6131,      6132,      6133,      6134,      6135,
      6136,      6137,      6138,      6139,      6140,      6141,      6142,      6143,      6144,      6145,      6146,      6147,      6148,      6149,      6150,      6151,
      6152,      6153,      6154,      6155,      6156,      6157,      6158,      6159,      6160,      6161,      6162,      6163,      6164,      6165,      6166,      6167,
      6168,      6169,      6170,      6171,      6172,      6173,      6174,      6175,      6176,      6177,      6178,      6179,      6180,      6181,      6182,      6183,
      6184,      6185,      6186,      6187,      6188,      6189,      6190,      6191,      6192,      6193,      6194,      6195,      6196,      6197,      6198,      6199,
      6200,      6201,      6202,      6203,      6204,      6205,      6206,      6207,      6208,      6209,      6210,      6211,      6212,      6213,      6214,      6215,
      6216,      6217,      6218,      6219,      6220,      6221,      6222,      6223,      6224,      6225,      6226,      6227,      6228,      6229]

    
############################################################################################################################
############################################################################################################################
i=1
for items in dof_map:
    node_number= items[1]
    











































    