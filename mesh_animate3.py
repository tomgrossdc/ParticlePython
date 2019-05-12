#mesh_animate.py



import numpy as np
from math import cos, sin, atan2, sqrt
import datetime, os, sys, time

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class mesh_animate(object):
    def __init__(self,lonlatbox=(-77.,36.,-75.,40.0)):

        # Create new Figure and an Axes which fills it.
        figwidth=8.
        yperx=(lonlatbox[3]-lonlatbox[1])/(lonlatbox[2]-lonlatbox[0])
        meterperlat = 111132.92 -559.82*cos(2.*lonlatbox[1])+1.175*cos(4.*lonlatbox[1])
        meterperlon = 111132.92
        figtall = figwidth*yperx*meterperlat/meterperlon
        # correct figsize for longitude latitude ratio
        
        self.fig = plt.figure(figsize=(figwidth,figtall))
        #self.ax = self.fig.add_axes([0,0,1,1], frameon=True)
        self.ax = self.fig.add_axes([.1,.1,.8,.8], frameon=True)
        self.lonlatbox=lonlatbox
        self.ax.set_xlim(self.lonlatbox[0],self.lonlatbox[2])
        #self.ax.set_xticks([])
        self.ax.set_ylim(self.lonlatbox[1],self.lonlatbox[3])
        #self.ax.set_yticks([])

        self.n_particles = 20
        self.n_time = 250
        #self.get_particles_time()

    def get_particles_pp(self,pp,pcolors,timep):
#         def get_particles_pp(self,pp,pcolors):
        # p is particle array from main_mesh.py
        # pcolors is output from MG.plot_box_pick
        self.n_particles = len(pcolors)
        self.n_time = np.shape(pp)[0]
        print("get_particles_pp, n_particles={0}, n_time={1}, pp={2}".format(self.n_particles, self.n_time,np.shape(pp)))
        self.particles=pp
        self.timep=timep

        colorlist=np.array(((1,0,0,1),(0,1,0,1),(0,0,1,1),(0,1,1,1),
                            (.5,.5,1,1),(0,.5,1,1),(0.5,0,.5,1),(0,0,0,1)))
        lenclist=len(colorlist)
        self.colorspp=np.zeros((self.n_particles,4))
        #print("pcolors={0}".format(pcolors))
        for i in range(self.n_particles):
            self.colorspp[i]=colorlist[pcolors[i]%lenclist]

        #self.colorspp=(1,0,0,1)
        
        
    def get_particles_time(self):
        # Create three dimensional array particles[n_time,n_particles,2]
        self.particles = np.zeros((self.n_time,self.n_particles,2))
        colorlist=np.array(((1,0,0,1),(0,1,0,1),(0,0,1,1),(0,1,1,1)))
        # Initialize the particles in box positions and with
        for frame_number in range(self.n_time):
            facto= (float(frame_number)%float(self.n_time))/float(self.n_time)
            xf=(0.1-.5)*facto + .5
            xl=(0.9-.5)*facto + .5
            yf=(0.1-.5)*facto + .5
            yl=(0.9-.5)*facto + .5
            x=np.linspace(xf,xl,self.n_particles/5)
            y=np.linspace(yf,yl,5)
            xv,yv=np.meshgrid(x,y,indexing='xy')
            xvv=np.reshape(xv,(1,self.n_particles))
            yvv=np.reshape(yv,(1,self.n_particles))
            zz=np.zeros((2,self.n_particles))
            zz[0]=xvv
            zz[1]=yvv
            zzz=np.transpose(zz)
        
            self.particles[frame_number] = zzz

            self.colorspp=np.zeros((self.n_particles,4))
            for ix in range(self.n_particles):
                if xvv[0,ix]>(0.75-.5)*facto + .5:  self.colorspp[ix] = colorlist[0]
                elif xvv[0,ix]>(0.5-.5)*facto + .5:  self.colorspp[ix] = colorlist[1]
                elif xvv[0,ix]>(0.25-.5)*facto + .5:  self.colorspp[ix] = colorlist[2]
                else :  self.colorspp[ix] = colorlist[3]
                


    def run_animation(self,xmesh=False,setinterval=10):

        if xmesh:
            #plt.triplot(xmesh.nodes[:,0], xmesh.nodes[:,1], xmesh.triwater, color=(.5,.5,1,.5))
            Coast=np.argwhere(xmesh.mask==10)
            plt.plot(xmesh.lon[Coast],xmesh.lat[Coast],'k.',markersize=3)
            self.ax.set_xlim(self.lonlatbox[0],self.lonlatbox[2])
            self.ax.set_ylim(self.lonlatbox[1],self.lonlatbox[3])
        else:
            self.ax.set_xlim(np.min(self.particles[:,:,0]),np.max(self.particles[:,:,0]))
            self.ax.set_ylim(np.min(self.particles[:,:,1]),np.max(self.particles[:,:,1]))
        
        plt.title("time today{0}".format((time.time())))
        # Construct the scatter which we will update during animation

        self.scat = self.ax.scatter(self.particles[0][:, 0], self.particles[0][:, 1],
                  s=10, lw=0.5, facecolors=self.colorspp)
        #self.scat.set_edgecolors(self.colorspp)
        
        # Construct the animation, using the update function as the animation director.
        animation = FuncAnimation(self.fig, self.update, interval=setinterval)
        plt.show()


    def update(self,frame_number):
        # Get an index which we can use to access time slices
        current_index = frame_number % self.n_time

        # Update the scatter collection, with the new colors and positions.
        #self.scat.set_edgecolors(self.particles[current_index]['color'])
        self.scat.set_offsets(self.particles[current_index])
        datestart=datetime.datetime(2016,1,1,0,0,0)   # ROMS basedate
        datetimep=datestart+datetime.timedelta(days=self.timep[current_index])
        self.ax.set_title("date time {0} {1}".format(datetimep.date(), datetimep.time()))

if __name__ == '__main__':
# mesh_animate(n_particles=20,n_time=250)
    manim=mesh_animate(25,250)
    manim.get_particles_time()
    manim.run_animation(False,200)
