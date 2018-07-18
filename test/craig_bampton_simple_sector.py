import amfe
import os, sys
import numpy as np
from scipy.sparse import linalg as splinalg
import scipy.sparse as sparse

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


def craig_bampton(M, K, master_dofs, slave_dofs, no_of_modes=5):
    '''
    Computes the Craig-Bampton basis for the System M and K with the input
    Matrix b.

    Parameters
    ----------
    M : ndarray
        Mass matrix of the system.
    K : ndarray
        Stiffness matrix of the system.
     master_dofs : ndarray
        input with dofs of master nodes
    slave_dofs : ndarray
        input with dofs of slave nodes
    no_of_modes : int, optional
        Number of internal vibration modes for the reduction of the system.
        Default is 5.
    one_basis : bool, optional
        Flag for setting, if one Craig-Bampton basis should be returned or if
        the static and the dynamic basis is chosen separately

    Returns
    -------
    if `one_basis=True` is chosen:

    V : array
        Basis consisting of static displacement modes and internal vibration
        modes

    if `one_basis=False` is chosen:

    V_static : ndarray
        Static displacement modes corresponding to the input vectors b with
        V_static[:,i] being the corresponding static displacement vector to
        b[:,i].
    V_dynamic : ndarray
        Internal vibration modes with the boundaries fixed.
    omega : ndarray
        eigenfrequencies of the internal vibration modes.

    Examples
    --------
    TODO

    Notes
    -----
    There is a filter-out command to remove the interface eigenvalues of the
    system.

    References
    ----------
    TODO

    '''
    # boundaries
    ndof = M.shape[0]       
    K_tmp = K.copy()
    
    K_bb = K_tmp[np.ix_(master_dofs, master_dofs)]
    K_ii = K_tmp[np.ix_(slave_dofs,  slave_dofs)]
    K_ib = K_tmp[np.ix_(slave_dofs,  master_dofs)]
    Phi = splinalg.spsolve(K_ii,K_ib)

    # inner modes
    M_tmp = M.copy()
    # Attention: introducing eigenvalues of magnitude 1 into the system
    M_bb = M_tmp[np.ix_(master_dofs, master_dofs)]
    M_ii = M_tmp[np.ix_(slave_dofs,  slave_dofs)]
    M_ib = M_tmp[np.ix_(slave_dofs,  master_dofs)]
    

    omega, V_dynamic = splinalg.eigsh(K_ii, no_of_modes, M_ii)


    num_of_masters = len(master_dofs)
    num_of_slaves = ndof - num_of_masters

    I = np.identity(num_of_masters)
    Zeros = np.zeros( (num_of_masters, no_of_modes))

    T_local_row_1 = np.hstack((I,Zeros))
    T_local_row_2 = np.hstack((Phi.todense(),V_dynamic))
    T_local = np.vstack((T_local_row_1,T_local_row_2))
    
    local_indexes = []
    local_indexes.extend(master_dofs)
    local_indexes.extend(slave_dofs)
    #redution_index = []
    #redution_index.extend(master_dofs)
    #redution_index.extend(list(range(num_of_masters+1,+num_of_masters+no_of_modes)))

    T = np.zeros((ndof, num_of_masters + no_of_modes))
    for i_local, i_global in enumerate(local_indexes):
       T[i_global,:] = T_local[i_local,:]




    #omega = np.sqrt(omega)
    
    return T



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

master_dofs, slave_dofs = boundary_dofs, interior_dofs

K_global = read_calculix.read_matrix(global_K_file)
M_global = read_calculix.read_matrix(global_mass_file)
Kii_stress = read_calculix.read_matrix(interface_K_stress_file)

K_stress = assemble_K_stress(K_global,Kii_stress,local2global_dict)
K_stress_sym = 0.5*(K_stress + K_stress.T) 
T = craig_bampton(M_global, K_stress_sym, master_dofs, slave_dofs , no_of_modes=5)

teste= 1

K_cg = T.T.dot(K_stress_sym.todense()).dot(T)
M_cg = T.T.dot(M_global.todense()).dot(T)
teste = 1