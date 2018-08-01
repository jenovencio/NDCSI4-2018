# importing lib and setting a list of meshes to be tested
import sys 
import amfe
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import linalg as splinalg
import copy
import mpl_toolkits.mplot3d as a3

# making the path relative to the project
local_folder = os.getcwd()[:]
os.chdir('..')
project_folder = os.getcwd()[:]
os.chdir(local_folder)
sys.path.append(r'project_folder')

mesh_inp = os.path.join(project_folder,r'data\ansys_simple_blade_disc_files\simple_blade_disc.inp')


m = amfe.Mesh()
m.import_inp(mesh_inp,1.0)



my_comp = amfe.CraigBamptonComponent()
my_comp.set_mesh_obj(m)

my_material = amfe.KirchhoffMaterial(E=210E9, nu=0.3, rho=7.86E3, plane_stress=True, thickness=1.0)
my_comp.set_domain('SOLID_1_1_SOLID_ELSET',my_material)


#-----------------------------------------------------------------------------------------------------
if False:
    K, f = my_comp.assembly_class.assemble_k_and_f()
    K_, f = my_comp.assembly_class.assemble_k_and_f_neumann()
    M = my_comp.assembly_class.assemble_m()
    amfe.save_object(K,'K.pkl')
    amfe.save_object(M,'M.pkl')
    amfe.save_object(f, 'f.pkl')
    amfe.save_object(my_comp,'my_comp.pkl')
else:
    K = amfe.load_obj('K.pkl')
    M = amfe.load_obj('M.pkl')
    f = amfe.load_obj('f.pkl')
    my_comp = amfe.load_obj('my_comp.pkl')


dirsub = m.get_submesh('phys_group', 'DIRICHLET_ELSET')
cyclic_low = m.get_submesh('phys_group', 'LOW_ELSET')
cyclic_high = m.get_submesh('phys_group', 'HIGH_ELSET')
cyclic_top_low = m.get_submesh('phys_group', 'TOP_LOW_ELSET')
cyclic_top_high = m.get_submesh('phys_group', 'TOP_HIGH_ELSET')



def plot_normal(submesh, ax=None, scale = 1000 ):

    if ax==None:
        ax = a3.Axes3D(plt.figure()) 

    node_list = submesh.create_node_list()
    X = []
    Y = []
    Z = []
    U = []
    V = []
    W = []
    for node_id in node_list:
        node_normal_vec = submesh.get_normal_to_node(node_id, method = 'average')
        node_coord = submesh.get_node_coord(node_id)
        (x,y,z) = node_coord*scale
        (dx,dy,dz) = node_normal_vec*scale
        u,v,w = (x + dx, y + dy, z + dz)
        X.append(x)
        Y.append(y)
        Z.append(z)
        U.append(u)
        V.append(v)
        W.append(w)
    
    ax.quiver(X, Y, Z, U, V, W, length=5, normalize = True, color = 'black')
    return ax



virtual_cyclic_top_high = copy.deepcopy(cyclic_top_high)
virtual_cyclic_top_high.rot_z(45)

def find_node_pairs(cyclic_top_low,virtual_cyclic_top_high, tol_radius = 1e-6):
    ''' find node pairs for contact given two submeshs

    parameters:
        cyclic_top_low : SubMesh
            SubMesh with the Master nodes 

        virtual_cyclic_top_high: SubMesh
            SubMesh with the Slaves nodes 

        tol_radius : float
            tolerance for finding node pairs, if a node pair do not respect the minimum 
            tolerance it will not considered as node pairs

        return : 
            contact_elem_dict : dict
                dict that poitns master nodes to slaves

    '''

    master_nodes = cyclic_top_low.create_node_list()
    slaves_nodes = virtual_cyclic_top_high.create_node_list()
    
    # master points to slave # master is a key and slave is value
    contact_elem_dict = {}
    for master_node in master_nodes:
        master_coord = cyclic_top_low.get_node_coord( master_node)
        min_dist = 1E8
        for slave_node in slaves_nodes:
            slave_coord = virtual_cyclic_top_high.get_node_coord(slave_node)
            dist = np.linalg.norm(master_coord - slave_coord)
            if dist<min_dist:
                slave_pair = slave_node
                min_dist = dist

        if min_dist>tol_radius:
            print('It was not possible to find a slave node for master node %i ')
        else:
            contact_elem_dict[master_node] = slave_node
    
    return contact_elem_dict


contact_elem_dict = find_node_pairs(cyclic_top_low,virtual_cyclic_top_high, tol_radius = 1e-6)

ax1 = a3.Axes3D(plt.figure()) 
amfe.plot3D_submesh(cyclic_top_high,ax=ax1, alpha=0.5, color='grey', plot_nodes=True, interface_nodes=True, scale = 1000)
amfe.plot3D_submesh(cyclic_top_low,ax=ax1, alpha=0.5, color='blue', plot_nodes=True, interface_nodes=True, scale = 1000)
amfe.plot3D_submesh(virtual_cyclic_top_high,ax=ax1, alpha=0.5, color='orange', plot_nodes=True, interface_nodes=True, scale = 1000)

ax1 = plot_normal(cyclic_top_high, ax=ax1, scale = 1000 )
ax1 = plot_normal(cyclic_top_low, ax=ax1, scale = 1000 )
plt.show()