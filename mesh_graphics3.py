"""
mesh_graphics.py
Code to generate a 2d cloud of random points, triangulate it,
# assign velocities to the nodes, interpolate within triangles
# and return velocities for a group of drogues.

# Add Roms netcdf reading functions.  Will refactor to separate files soon.

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from math import cos, sin, atan2, sqrt

import os

DEBUG=False

def plot_box_pick(nx,ny,lonlatbox):
    # make an array of p at lower corner xo,yo, upper right x1,y1
    # lonlatbox=(x0=-76.18,y0=37.7,x1=-76.,y1=37.73)
    x0=lonlatbox[0]
    y0=lonlatbox[1]
    x1=lonlatbox[2]
    y1=lonlatbox[3]
    # pcolors is now just an index into several groups 1:6 or so
    pcolors=np.zeros(nx*ny,'int')
    p=np.zeros((nx*ny,2))
    numparticles=nx*ny
    ii=0
    for ix in range(nx):
        for iy in range(ny):
            p[ii,0]=x0+ix*(x1-x0)/(nx-1.)
            p[ii,1]=y0+iy*(y1-y0)/(ny-1.)
            pcolors[ii]=((ix/5)%6 + (0*iy/5)%6)%6
            ii+=1

    return p,pcolors,numparticles        

def plot_mesh_pick_line(xmesh):
    # in mainline do: 
    #plot_mesh2(xmesh)
    #
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #plot triangle meshes and particle initial points
    #ax.triplot(xmesh.nodes[:,0], xmesh.nodes[:,1], xmesh.tri.simplices, color='green')
    ax.plot(xmesh.nodes[:,0], xmesh.nodes[:,1], '.', color='lightblue')
    line, = ax.plot([np.mean(xmesh.nodes[:,0])],[np.mean(xmesh.nodes[:,1])])
    linebuilder = PointBuilder(line)
    plt.show()
    print ("x,y",linebuilder.xs[-1],linebuilder.ys[-1])
    p = np.transpose(np.array([linebuilder.xs[1:],linebuilder.ys[1:]]))
    return p

def plot_mesh_pick_line2(xmesh,numtimes=12):
    # in mainline do: 
    #p=plot_mesh_pick_line2(xmesh,numtimes=12)
    # p[numparticles,numtimes,2]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #plot triangle meshes and particle initial points
    #ax.triplot(xmesh.nodes[:,0], xmesh.nodes[:,1], xmesh.tri.simplices, color='green')
    ax.plot(xmesh.nodes[:,0], xmesh.nodes[:,1], '.', color='lightblue')
    line, = ax.plot([np.mean(xmesh.nodes[:,0])],[np.mean(xmesh.nodes[:,1])])
    linebuilder = PointBuilder(line)
    plt.show()
    if DEBUG: print ("x,y",linebuilder.xs[-1],linebuilder.ys[-1])
    p = np.transpose(np.array([linebuilder.xs[1:],linebuilder.ys[1:]]))
# get rid of the time index for this old version of p
#    pp=np.zeros((np.shape(p)[0],numtimes,2))
#    for i in range(np.shape(p)[0]):
#        pp[i,0]=p[i]
    numparticles = np.shape(p)[0]
    
    #pcolors is just a 0:6 index for groups of particles used by plotting to pick colors later
    pcolors=np.zeros(numparticles,'int')
    for i in range(numparticles):
        pcolors[i]=i%6
    return p,pcolors,numparticles
    
class PointBuilder:
    def __init__(self,line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        #print('click', event)
        if event.inaxes!=self.line.axes: return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()
        #self.point(self.xs,self.ys,'.',color="yellow")
        if DEBUG: print("xs, ys",self.xs[-1],self.ys[-1])


def plot_initialize():
  #Set up display window for plotting
    import matplotlib as mpl
    mpl.rcParams['figure.figsize'] = [9.0, 9.0]
    mpl.rcParams['figure.dpi'] = 80
    mpl.rcParams['savefig.dpi'] = 100

    mpl.rcParams['font.size'] = 16
    mpl.rcParams['legend.fontsize'] = 'large'
    mpl.rcParams['figure.titlesize'] = 'medium'
    print ("start plot_initialize")
   
# clear out the first call to fig
    plt.close()




def plot_particles(xmesh,p,pcolors):
    plt.triplot(xmesh.nodes[:,0], xmesh.nodes[:,1], xmesh.triwater, color='lightblue')
    colordays=['red','blue','green','orange','pink','purple']
    lenp=np.shape(p)[1]
    for ip in range(np.shape(p)[1]):
        plt.plot(p[:,ip,0],p[:,ip,1],color=colordays[min(5,pcolors[ip])],linestyle='-',linewidth=1)
        #plt.scatter(p[:,ip,0],p[:,ip,1],c=colordays[])
    plt.show()

def plot_particles_spots(xmesh,p,pcolors):
    plt.triplot(xmesh.nodes[:,0], xmesh.nodes[:,1], xmesh.triwater, color='lightblue')
    colordays=['red','blue','green','orange','pink','purple']
    lenp=np.shape(p)[1]
    istepper=12
    for it in range(0,np.shape(p)[0],istepper):
        plt.plot(p[it,:,0],p[it,:,1],'o',color=colordays[int(it/istepper)%len(colordays)])
        #plt.scatter(p[:,ip,0],p[:,ip,1],c=colordays[])
    plt.show()



def plot_mesh(xmesh, ymesh):

    #fig, ax = plt.subplots()

    #plot triangle meshes 
    plt.triplot(xmesh.nodes[:,0], xmesh.nodes[:,1], xmesh.tri.simplices, color='brown')
    #plt.triplot(xmesh.nodes[:,0], xmesh.nodes[:,1], xmesh.triwater, color='lightblue')
    #plt.plot(xmesh.nodes[:,0], xmesh.nodes[:,1], '.', color='red')

    #plt.triplot(ymesh.nodes[:,0], ymesh.nodes[:,1], ymesh.tri.simplices, color='yellow')
    plt.triplot(ymesh.nodes[:,0], ymesh.nodes[:,1], ymesh.triwater, color='blue')
    #plt.plot(ymesh.nodes[:,0], ymesh.nodes[:,1], '.', color='orange')

    #plt.plot(p[:,0],p[:,1],'+')
    #cid = fig.canvas.mpl_connect('button_press_event', onpick) 
    #fig.canvas.mpl_disconnect(cid)

    plt.show()
    #print (" x=%d , y=%d ", x_mouse,y_mouse)
    #print ("p", p.shape)
    #print(p)
    
    #pelements = xmesh.tri.find_simplex(p[:,:])
    #print ("pelements ")
    #print (pelements)
    #puniq= np.unique(pelements)
    #print puniq




def plot_velocities ( xmesh, U1,V1, p,Up,Vp):
    # plot_velocities(xmesh,U24[4],U24[4],p,Up,Vp)
    DEBUG = True
    #mesh plot of velocities
    if DEBUG :
        plt.figure()
        plt.gca().set_aspect('equal')
        plt.tricontourf(xmesh.nodes[:,0],xmesh.nodes[:,1],xmesh.tri.simplices,U1)
        plt.colorbar()
        plt.title('Contour plot of U1 on xmesh')
        #plt.quiver(p[:,0],p[:,1],Up,Vp,units='width',width=0.003,pivot='tail',color='red')
        plt.show()
        plt.figure()
        plt.gca().set_aspect('equal')
        plt.tricontourf(xmesh.nodes[:,0],xmesh.nodes[:,1],xmesh.tri.simplices,V1)
        plt.colorbar()
        plt.title('Contour plot of V1 on xmesh')
        plt.show()


    """
    Q2=plt.quiver(xmesh.nodes[:,0], xmesh.nodes[:,1],Umesh,Vmesh,units='width',width=0.003,pivot='tail',color='red')
    plt.show
    """


    
"""
    # plot velocity vectors at particle points and ymesh triangulation
    Q2=plt.quiver(p[:,0],p[:,1],Urot,Vrot,units='width',width=0.003,pivot='tail',color='black')
    plt.triplot(ymesh.nodes[:,0], ymesh.nodes[:,1], ymesh.tri.simplices, color='yellow')
    plt.show()
    #Q2=plt.quiver(p[:,0],p[:,1],Urot,Vrot,units='width',width=0.003,pivot='tail',color='black')
    #plt.triplot(ymesh.nodes[:,0], ymesh.nodes[:,1], ymesh.tri.simplices, color='orange')
    #plt.show()
"""



