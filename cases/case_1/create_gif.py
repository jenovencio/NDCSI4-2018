# importing lib and setting a list of meshes to be tested
import sys 
import amfe
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import linalg as splinalg
import copy
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import mpl_toolkits.mplot3d as a3

# making the path relative to the project
local_folder = os.getcwd()[:]
os.chdir('..\..')
project_folder = os.getcwd()[:]
os.chdir(local_folder)
sys.path.append(project_folder)


from utils.write import *


try:

    m = amfe.load_obj('mesh.pkl')
except:
    mesh_inp = os.path.join(project_folder,r'cases\case_1\simple_blade_disc.inp')
    m = amfe.Mesh()
    m.import_inp(mesh_inp,1.0)
    amfe.save_object(m,'mesh.pkl')

angle = 45
meshes = [m]
for a in np.arange(angle,360,angle):
    meshes.append(m.rot_z(a))


def update(i):
    for m in meshes:
        ax = amfe.plot3Dmesh(m, ax=ax1, boundaries=True, alpha=0.2, color='blue', scale = 1000, plot_nodes=True)
    
    ax.view_init(0 + i, 0 )
    plt.draw()
    ax.set_aspect('equal', adjustable='box')
    #xlim = ax.get_xlim()
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    zlim = ax.get_zlim()

    min_ = min(xlim[0],ylim[0],zlim[0])
    max_ = max(xlim[1],ylim[1],zlim[1])
    lim = (min_,max_)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_zlim(lim)

    #ax.set_zlim(xlim)
    plt.axis('off')
    return ax

def create_gif():
    # FuncAnimation will call the 'update' function for each frame; here
    # animating over 10 frames, with an interval of 200ms between frames.
    Writer = animation.writers['imagemagick']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

    anim = FuncAnimation(fig, update, frames=np.arange(0, 180, 5), interval=200)
    anim.save('geometry.gif', dpi=400, writer=writer)


if __name__ == '__main__':
    fig = plt.figure()
    ax1 = a3.Axes3D(fig) 
    #ax1.view_init(90, -90)
    create_gif()
    
