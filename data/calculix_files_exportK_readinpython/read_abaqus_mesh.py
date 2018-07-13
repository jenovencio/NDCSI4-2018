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
    
    return 
    
    nodes_dict: dict
        dict which maps nodes to coordinates
    elem_list: list
        list with elements sets and element to node mapping 
        
    nset_list : list
        list with node sets
    
    '''
    
    file_labels = {}    
    file_key_words = ['ELEMENT','NSET','NODE']
    with open(filepath) as f:
        for line in f:            
            try:
                words_list = line.split(',')
                first_word = words_list[0].replace(',','').replace('\n','')
            except:
                continue
            

            if first_word[1:] == file_key_words[0]:
                elem_type = words_list[1].split('=')[1].replace(' ','')
                elem_set = words_list[2].split('=')[1].replace(' ','').replace('\n','')
                dict_key = (first_word[1:], elem_type, elem_set) 
                
            elif first_word[1:] == file_key_words[1]: 
                nset_key = words_list[1].split('=')[1].replace(' ','').replace('\n','')
                dict_key = (first_word[1:], nset_key) 
            
            elif first_word[1:] == file_key_words[2]:
                dict_key = first_word[1:] 
            
            try:
                if dict_key not in file_labels:
                    file_labels[dict_key] = []   
            except:
                continue
                    
            file_labels[dict_key].append([line])

            
    # parsing list to to dict date structure   
    elem_list = []
    nset_list = []    
    for key,values in file_labels.iteritems():
        if file_key_words[0] in key:
            
            elem_set_data = {}
            elem_set_data['elem_type'] = key[1]
            elem_set_data['elem_set'] = key[2]
            elem_set_data['elem_dict'] = parse_elem_list(values)
            elem_list.append(elem_set_data)
            
        elif file_key_words[1] in key:
            elem_set_data = {}
            elem_set_data['node_set'] = key[1]
            elem_set_data['node_list'] = parse_nset(values)
            nset_list.append(elem_set_data)
            
        elif file_key_words[2] in key:
            nodes_dict = parse_list(values)
            
        else:
            raise('Error parsing file')
    
    
    return nodes_dict, elem_list, nset_list 
    
    
def parse_list(list_with_strings):
    ''' receives a string list with nodes and
        return a dict with maps node number to coordinates
    
    argument:
       list_with_strings : list
    
    return:
        local_dict : dict
            dict with node as keys as coordinates as values
    '''
    
    local_dict = {}
    for line in list_with_strings[1:]:
        values = line[0].split(',')        
        local_dict[int(values[0])] = [float(i) for i in values[1:]]
    return local_dict     
    
def parse_elem_list(list_with_strings):
    ''' receives a string list with elem set and nodes information
    and creates a dict mapping element to nodes 
    
    argument:
        list_with_strings : list
    
    return:
        local_dict : dict
            dict with elements as keys as nodes as values
        
    
    '''
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
    ''' receives a string list of node set information
    and creates a node list 
    
    argument:
        nsets_str_list : list
    
    return:
        node_list : list
        
    
    '''
    for line in nsets_str_list[1:]:
        try:
            values = line[0].replace('\n','').split(',')    
            if line[0][-2] == ',':
                node_list = [int(i) for i in values[1:-1]]           
            else:  
                node_list.extend([int(i) for i in values])
                
        except:
            continue
        
    return node_list
    
    

if __name__ == "__main__":
    import os
    local = os.chdir('..')
    filepath = r'data\beampl.inp'
    nodes_dict, elem_list, nset_list  = read_inp(filepath)
    

