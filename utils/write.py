def writeMatrix(M,filename):

    ndof, garbage = M.shape
    with open(filename,'w') as f:

        f.write(str(ndof)) 
        f.write('\n') 
        for columns in M.todense().T:
            for i in columns.tolist()[0]:
                f.write(str(i))
                f.write('\n') 
    return None
                
def write_map_matrix(filename,inv_id_matrix,local2global_master, label=True):

    with open(filename,'w') as f:
        
        f.write(str(2 + len(local2global_master)))
        f.write('\n')
        
        if label:
            f.write('Node dir dof')
            f.write('\n')
        
        for key, global_dof in local2global_master.items():
            global_node = inv_id_matrix[global_dof]['node']
            direction = inv_id_matrix[global_dof]['direction']
            f.write(str(global_node) + ' ' + str(direction) + ' ' + str(key))
            f.write('\n')
    
    return None
    

def write_contact(filename,contact_obj):
    ''' write setup file given a contact object
    
    '''
    with open(filename,'w') as f:
        num_of_contact_elem  = len(contact_obj.contact_elem_dict)
        num_of_lines = 2 + (2 + 3)*num_of_contact_elem
        
        # writing number of lines
        f.write(str(num_of_lines) + '\n')
        
        # writing number of contact elements
        f.write(str(num_of_contact_elem)+ '\n')
        
        # writing master and slaves
        for master_id, slave_id in contact_obj.contact_elem_dict.items():
            f.write(str(master_id) + '\n')
            f.write(str(slave_id) + '\n')
        
        for master, normal_vec in contact_obj.master_normal_dict.items():
            for component in normal_vec:
                f.write(str(component) + '\n')
    
    return None