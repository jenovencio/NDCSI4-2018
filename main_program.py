# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 15:51:29 2018

@author: SIV0021
"""

"""
read mesh file in .inp format
"""

mesh_file = open("beampl.inp","r")


file = open("testfile.txt","w") 
 
file.write("Hello World1123") 

file.close() 

fid = fopen([beampl,'.inp']);


CellRead = textscan(fid, '%s', 'Delimiter', '\n');


fclose(fid)
