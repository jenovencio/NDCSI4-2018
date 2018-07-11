# -*- coding: utf-8 -*-
"""


reading input file from Abaqus

Author: Guilherme Jenovencio, M.Sc.
email: guilherme.jenovencio@gmail.com
"""

#import pandas as pd


def read_inp(filepath):
    ''' This function parse Abaqus input file
    
    arguments
        filepath : str
            path of Abaqus filename
    
    return nodes and element dataframe
    '''
    
    file_labels = {}
    file_labels['*NODE'] = []
    file_labels['*ELEMENT'] = []
    file_labels['*NSET'] = [] 
    
    file_keys = file_labels.keys()
    count = 0
    with open(filepath) as f:
        for line in f:
            
            try:
                first_word = line.split()[0].replace(',','')
            except:
                continue
            
            if first_word in file_keys:
                print(first_word)
                local_list = file_labels[first_word]
                count+=1
            
            if count>0:   
                local_list.append([line])
            
    # parsing nodes to numpy array        
    nodes_dict = parse_list(file_labels['*NODE'])
    elem_dict = parse_elem_list(file_labels['*ELEMENT'])
    nset_dict = parse_nset(file_labels['*NSET'])
    
    return nodes_dict, elem_dict, nset_dict 
    
    
def parse_list(list_with_strings):
    local_dict = {}
    for line in list_with_strings[1:]:
        values = line[0].split(',')        
        local_dict[int(values[0])] = [float(i) for i in values[1:]]
    return local_dict     
    

def parse_elem_list(list_with_strings):
    local_dict = {}
    for line in list_with_strings[1:]:
        values = line[0].split(',')        
        if line[0][-2] == ',':
            node_list = [int(i) for i in values[1:-1]]
            elem_key  = int(values[0])
        else:  
            node_list.extend([int(i) for i in values])
            local_dict[elem_key] = node_list
    
    return local_dict       
    

def parse_nset(nsets_str_list):
    local_key = nsets_str_list[0][0].split(',')[1].split('=')[1].replace('\n','')
    local_dict = {}
    for line in nsets_str_list[1:]:
        values = line[0].split(',')        
        if line[0][-2] == ',':
            node_list = [int(i) for i in values[1:-1]]           
        else:  
            node_list.extend([int(i) for i in values])
            local_dict[local_key] = node_list
    
    return local_dict   
    
    

#if '__name__'=='__main__':
filepath = r'C:\Projects\data\beampl.inp'
nodes_dict, elem_dict, nset_dict  = read_inp(filepath)

