import amfe
import os, sys
import numpy as np
# making the path relative to the project
local_folder = os.getcwd()[:]
os.chdir('..')
project_folder = os.getcwd()[:]
os.chdir(local_folder)

wrappers_folder = os.path.join(project_folder , 'wrappers')
sys.path.append(wrappers_folder )

import read_abaqus_matrix as read_calculix


def find_interior_boundary_dofs(global_id_matrix,local_id_matrix):
    ''' convert local to global

    '''
    interior_dofs = []
    all_dofs = []
    boundary_dofs = []
    local2global_dict = {}
    global2local_dict = {}

    for global_node in global_id_matrix:
        all_dofs.extend(global_id_matrix[global_node])
        if global_node in local_id_matrix:
            if len(global_id_matrix[global_node])==len(local_id_matrix[global_node]):
                ''' because there no information of dof and x,y and z coord 
                '''
                interior_dofs.extend(global_id_matrix[global_node])
                for i,dof in enumerate(global_id_matrix[global_node]):
                    local_dof = local_id_matrix[global_node][i]
                    local2global_dict[local_dof] = dof
                    global2local_dict[dof] = local_dof

            else:
                print('ERROR')
        else: 
            boundary_dofs.extend(global_id_matrix[global_node])

            

    #boundary_dofs = list(set(all_dofs) - set(interior_dofs))
    return interior_dofs, boundary_dofs, local2global_dict, global2local_dict




def assemble_K_stress(K_global,Kii_stress,local2global_dict):
    ''' 
        replace K_global by stiffness matrix

        arguments :
            K_global : scipy.sparse.matrix
                sparse matrix without stress
            Kii_stress : scipy.sparse.matrix
                sparse matrix without stress
            local2global_dict : dict
                dict which maps local dofs to global dofs

        returns 
            K_global_stress : scipy.sparse.matrix
                return a sparse matrix with stress effect
    '''
    for local_dof, global_dof in local2global_dict.items():
        K_global[global_dof,global_dof] = Kii_stress[local_dof,local_dof]
    
    return K_global







# files to build K and M matrix
data_folder = os.path.join(project_folder,r'data\calculix_simple_sector')

global_id_file = os.path.join(data_folder, 'case_global.dof')
global_K_file = os.path.join(data_folder,'case_global.sti')
global_mass_file = os.path.join(data_folder,'case_global.mas')
interface_K_stress_file = os.path.join(data_folder,'case_stress.sti')
interior_id_file = os.path.join(data_folder,'case_stress.dof')

global_id_matrix = read_calculix.read_dof_file(global_id_file)
local_id_matrix = read_calculix.read_dof_file(interior_id_file)

interior_dofs, boundary_dofs, local2global_dict, global2local_dict = find_interior_boundary_dofs(global_id_matrix,local_id_matrix)


K_global = read_calculix.read_matrix(global_K_file)
M_global = read_calculix.read_matrix(global_mass_file)
Kii_stress = read_calculix.read_matrix(interface_K_stress_file)

K_stress = assemble_K_stress(K_global,Kii_stress,local2global_dict)
K_stress_sym = 0.5*(K_stress + K_stress.T) 
T = amfe.craig_bampton(M_global, K_stress_sym, np.array(boundary_dofs), no_of_modes=5, one_basis=True)

T.shape