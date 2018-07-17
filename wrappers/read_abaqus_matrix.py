import sys
import numpy as np
import os
import scipy.sparse as sparse

# making the path relative to the project
local_folder = os.getcwd()[:]
os.chdir('..')
project_folder = os.getcwd()[:]
os.chdir(local_folder)


def read_dof_file(filename):
    ''' create a dictionary to map nodes to dofs
    based on filename

    arguments :
        filename : str
            path of the file containing dofs

    return:
        dict
            maps among nodes and dofs
            dict has a python index so 1 is 0
    '''
    ID_matrix = {}
    with open(filename,'r') as f:
        for i,line in enumerate(f):
            value = line.split(".");
            key = int(value[0]) - 1
            dof = i
            if key in ID_matrix: 
                ID_matrix[key].extend([dof])
            else:
                ID_matrix[key] = [dof]

    return ID_matrix

 
def read_matrix(filename,number_of_lines=None):
    ''' return a matrix based on a file

    arguments:
        filename : str
            path of the filename
        
       number_of_lines : int
           number of lines of the file
    
    returns
        A : np.array
            matrix based on the file with python indexes
    '''

    if number_of_lines is None:
        with open(filename,'r') as f:
            number_of_lines = len(f.readlines())

    row_list = []
    col_list = []
    data = []
    with open(filename,'r') as f:
        for line in f:
            values = line.split()
            i = int(values[0]) - 1
            j = int(values[1]) - 1
            value = float(values[2])
            row_list.append(i)
            col_list.append(j)
            data.append(value)

    A = sparse.csc_matrix((data,(row_list,col_list)))
    return A



if "__name__" == "__main__":
    """ Read the dof array  """
    dof_file = os.path.join(project_folder, r'data\files\id.dof')
    ID_matrix = read_dof_file(dof_file)

    """ read the mass matrix """
    
    mass_file =   os.path.join(project_folder,r'data\files\mass.mas')   
    M = read_matrix(mass_file)   


    """ read the stiffness matrix """
    sti_file = os.path.join(project_folder,r'data\files\k.sti')  

    K = read_matrix(sti_file)   

