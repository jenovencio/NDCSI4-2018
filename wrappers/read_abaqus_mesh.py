# -*- coding: utf-8 -*-
"""


reading input file from Abaqus

Author: Guilherme Jenovencio, M.Sc.
email: guilherme.jenovencio@gmail.com
"""

import pandas as pd
import numpy as np



# Same for Abaqus
abaq2amfe = {'C3D10M' : 'Tet10',
             'C3D8' : 'Hexa8',
             'C3D20' : 'Hexa20',
             'C3D20R' : 'Hexa20',
             'C3D4' : 'Tet4',
             'C3D6' : 'Prism6', # 6 node prism
             'C3D8I' : 'Hexa8', # acutally the better version
             'B31' : None,
             'CONN3D2' : None,
            }

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
    for key,values in file_labels.items():
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
    


def create_amfe_elem_data_frame(elem_list):
    '''
    '''
    columns_list=['idx_abaqus',
             'el_type',
             'no_of_tags',
             'phys_group',
             'geom_entity',
             'no_of_mesh_partitions',
             'partition_id',
             'partitions_neighbors']
    
    node_idx = len(columns_list)
    el_df = pd.DataFrame(columns=columns_list)
    
    for elem_set_tags in elem_list:
        elem_set = elem_set_tags['elem_set'] 
        elem_type = elem_set_tags['elem_type'] 
        my_dict = {}
        for elem in elem_set_tags['elem_dict']:             
            my_dict['idx_abaqus'] = elem
            my_dict['phys_group'] = elem_set
            my_dict['el_type'] = elem_type
            my_dict['no_of_mesh_partitions'] = 1
            my_dict['partition_id'] = 1
            node_list = elem_set_tags['elem_dict'][elem]
            for id_col, value in enumerate(node_list):
                # changing node index to amfe index starting from 0
                my_dict[str(id_col)] = value - 1
                if str(id_col) not in el_df:
                    el_df[str(id_col)] = 0
                
            el_df = el_df.append(pd.Series(my_dict),ignore_index=True)    
    
    # change the el_type to the amfe convention
    el_df['el_type'] = el_df.el_type.map(abaq2amfe)
    return el_df, node_idx 

def create_amfe_node_array(nodes_dict):   
    ''' convert dict to amfe nodes format
    
    argumets:
        nodes_dict : dict
            dict with number node to coordinates
    returns
        np.array
    '''
    
    node_list = np.sort(np.array(list(nodes_dict.keys())),axis=-1)
    node_coord = []
    for node in node_list:
        node_coord.append(nodes_dict[node])
        
    return np.array(node_coord)
        
if __name__ == "__main__":
    import os
    local = os.chdir('..')
    filepath = r'data\beampl.inp'
    nodes_dict, elem_list, nset_list  = read_inp(filepath)
    el_df, node_idx = create_amfe_elem_data_frame(elem_list)
    nodes = create_amfe_node_array(nodes_dict)

    import sys 
    sys.path.append(r'C:\AMfe')
    
    import amfe
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d as a3
    m2 = amfe.Mesh()
    m2.el_df = el_df
    m2.node_idx = node_idx
    m2.nodes = nodes
    
    fig = plt.figure(figsize=(30, 30), dpi= 30, facecolor='w', edgecolor='k')
    ax = a3.Axes3D(fig)
    ax = amfe.plot3Dmesh(m2,ax,alpha=0.2, plot_nodes=False)
    #ax.set_axis_off()
    #ax.view_init(90, -90)
    #ax.set_xlim([-80,80])
    #ax.set_ylim([90,250])
    #ax.set_zlim([-80,80])
    plt.show()
