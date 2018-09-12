# importing lib and setting a list of meshes to be tested
import sys 
import amfe
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import linalg as splinalg
import copy



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

virtual_cyclic_top_high = copy.deepcopy(cyclic_top_high)
virtual_cyclic_top_high.rot_z(45)

contact =amfe.Contact(cyclic_top_low, virtual_cyclic_top_high)

contact2 = amfe.Cyclic_Contact(cyclic_top_low, cyclic_top_high, sector_angle = 45)

