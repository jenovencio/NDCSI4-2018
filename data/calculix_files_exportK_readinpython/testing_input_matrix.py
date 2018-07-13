# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 12:18:34 2018

@author: SIV0021
"""

import sys
sys.path.append(r'D:\NDCSI\framework')

from wrappers import read_abaqus_mesh  as w

w.read_inp('MESH.inp')