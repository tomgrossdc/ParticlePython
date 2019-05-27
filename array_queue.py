#array_queue.py
# input the ncfilename and return the xmesh,ymesh,amesh data structures
# runs all three in parallel to speed up the startup of the ParticleTracker
# 90% of time is spent in the inefficient mask generator
#     xmesh,ymesh,amesh=array_queue(ncfilename)

# basic outline of particle tracker with multiprocessing

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.colors import ListedColormap
from math import cos, sin, atan2, sqrt


from netCDF4 import Dataset
from scipy.spatial import Delaunay
import time

from multiprocessing import Lock, Process, Queue, current_process
import queue
import multiprocessing as mp


# create mesh objects with self.lon, lat, mask
class xyarray(object):
    def __init__(self,lon,lat,mask,DataType,ModelType='ROMS_REGULAR'):
        # or ROMS_FIELDS  DataType= U V W ANGLE
        self.lon = lon
        self.lat = lat
        self.mask = mask
        self.DataType = DataType
        self.ModelType=ModelType
        self.s_rho = 0

    def reshape_mesh(self):
        self.dx=self.lon[2][3]-self.lon[2][2]
        self.dy=self.lat[3][2]-self.lat[2][2]
        #for i in range(25000000):   # one second for 25 million multiplications / processor
        #    x=i*i
        print(" lon, lat shapes ",self.lon.shape,self.lat.shape, self.DataType)
        print(" dx,dy ",self.dx,self.dy)
        #time.sleep(5)

# finds coast and reshapes to 1D array
        print("old mask reshape ",self.mask.shape[0]*self.mask.shape[1] )

        #numpy gradient method for mask on coast
        if self.ModelType == 'ROMS_REGULAR':
            self.mask_fix2()

        # Inefficient i,j differences
        #self.mask_fix()

        # No coast locations at all
        if self.ModelType == 'ROMS_FIELDS':
            self.mask = self.mask.reshape(self.mask.shape[0]*self.mask.shape[1])
            #self.mask=np.ones(self.mask.shape)
            
        self.lon = self.lon.reshape(self.lon.shape[0]*self.lon.shape[1])
        self.lat = self.lat.reshape(self.lat.shape[0]*self.lat.shape[1])
        print("shapes ",self.lat.shape,"mask",np.shape(self.mask))
        self.x_ = self.lon[self.mask>0]
        self.y_ = self.lat[self.mask>0]
        self.angle = 0.0 * self.lon[self.mask>0]
        self.nodes = np.transpose(np.array([self.x_, self.y_]))
        self.range = .2
        print("end of reshape_mesh",np.shape(self.nodes))        
    
    def mask_fix2(self):
        print(" Start mask_fix2 shape", self.mask.shape,self.DataType,self.ModelType)
        timefirst=time.time()
        maskfix=self.mask

        gradmask=np.gradient(self.mask)
        absgradmask=np.abs(gradmask)
        abstotalgrad=absgradmask[0]+absgradmask[1]
        N=self.mask.shape[0]*self.mask.shape[1]
        abstotalgrad=abstotalgrad.reshape(N)
        maskfix=maskfix.reshape(N)
        Coast=[]
        for ij in range(N):
            if maskfix[ij]==0 and abstotalgrad[ij]>0 :
                Coast.append(ij)
        maskfix[Coast]=10
        
        self.mask = maskfix
        print (" end mask_fix2",time.time()-timefirst,self.mask.shape )


    def mask_fix(self):
        print(" Start mask_fix shape", self.mask.shape,self.DataType,self.ModelType)
        timefirst=time.time()
        maskfix=np.zeros(self.mask.shape)
        # experiments with masking
        # sum nine mask points around a point
        # if ==0 and sum>0 it is a coast point set to 10.
        # else keep as is.
        #  massive problems with trying to equate a scratch array to self.mask
        # caused infinite sums and worse.
        # Also indexing was wierd.  [j-1:j+2] seems to sums j-1 to j+1 wierd
        #
        #iijj=0
        sumij=1.
        for i in range(1,(self.mask.shape[0]-1) ):
            #march across the columns in x,i direction
            #print(" 693 ",i,iijj,sumij)
            sumij0=0.
            sumij1=0.
            sumij2=0.
            for j in range(1,(self.mask.shape[1]-1) ):
                # crawl down the x(i) column from j=0:M  I don't get it
                sumij0=sumij1
                sumij1=sumij2
                sumij2=sum(self.mask[i+1][j-1:j+2])
                sumij = sumij0+sumij1+sumij2
                #sumij=sum(self.mask[i-1][j-1:j+2])+sum(self.mask[i][j-1:j+2])+sum(self.mask[i+1][j-1:j+2])
                if self.mask[i][j]<0.5 and sumij>1.2 :
                    # point is a masked value adjacent to a non-masked point
                    #iijj+=1
                    #print( iijj,i-1,i,i+1,j-1,j,j+1, sumij)
                    #print("   ",self.mask[i-1][j+1],self.mask[i][j+1],self.mask[i+1][j+1])
                    #print("   ",self.mask[i-1][j],self.mask[i][j],self.mask[i+1][j+1])
                    #print("   ",self.mask[i-1][j-1],self.mask[i][j-1],self.mask[i+1][j-1])
                    maskfix[i][j] = 10.
                elif self.mask[i][j]<0.5:
                    maskfix[i][j]=0.0
                else:
                    maskfix[i][j]=1.

        #self.maskfix = self.maskfix.reshape(self.mask.shape[0]*self.mask.shape[1])
        self.mask = maskfix.reshape(self.mask.shape[0]*self.mask.shape[1])
        # now we have a full maskfix turned 1D with 10 at points on coast
        # Should still have 1 at interior points
        print (" end mask_fix",time.time()-timefirst,self.mask.shape )




    def triangulate(self):
        #New Triangulate Fields
        #Read in file of outside of grid points
        #Remove all triangles containing those points from tri.simplices
        
        DEBUG=True
        print("triangulate ",self.ModelType,self.DataType)
    # call Delaunay triangulation to return
    # self.elements = [node1, node2, node3]
    # self.adjacent = [element1, element2, element3]
        self.tri = Delaunay(self.nodes)

    #def makefactors(self):
    # Calculate the linear interpolation factors
        if DEBUG : print ("makefactors",self.ModelType)
        NODES=self.tri.simplices
        ii=0
        icoast=0
        #self.triwater=self.tri.simplices
        triwatersn = np.arange(NODES.shape[0])
        ddiagnol = self.dx**2 + self.dy**2

        self.a=np.zeros(NODES.shape)
        self.b=np.zeros(NODES.shape)
        self.c=np.zeros(NODES.shape)
        for N in NODES:
            #factors for triangle ii
            X=self.nodes[N,0]
            Y=self.nodes[N,1]
            # WHY is this test inserted with "not",  temp change to disconnect
            #if  (self.ModelType=="ROMS_REGULAR"):
            if  True:
                d=max(X[0]*Y[1]+X[1]*Y[2]+X[2]*Y[0] - X[0]*Y[2]-X[1]*Y[0]-X[2]*Y[1], .000001)
                self.a[ii,0]=(Y[1]-Y[2])/d
                self.b[ii,0]=(X[2]-X[1])/d
                self.c[ii,0]=(X[1]*Y[2] - X[2]*Y[1])/d
                self.a[ii,1]=(Y[2]-Y[0])/d
                self.b[ii,1]=(X[0]-X[2])/d
                self.c[ii,1]=(X[2]*Y[0] - X[0]*Y[2])/d
                self.a[ii,2]=(Y[0]-Y[1])/d
                self.b[ii,2]=(X[1]-X[0])/d
                self.c[ii,2]=(X[0]*Y[1] - X[1]*Y[0])/d
            #print (ii,N,X,Y,d)
            # find only the small triangles
            d1=(X[0]-X[1])**2 + (Y[0]-Y[1])**2
            d2=(X[1]-X[2])**2 + (Y[1]-Y[2])**2
            d3=(X[0]-X[2])**2 + (Y[0]-Y[2])**2
            #print(" triangulate d1 1",sqrt(d1),sqrt(d2),sqrt(d3),sqrt(ddiagnol))
            if (d1+d2+d3)<(6.*ddiagnol) :
                #print(" %f+%f+%f less than %f"%(sqrt(d1),sqrt(d2),sqrt(d3),sqrt(ddiagnol)))
                #self.triwater[icoast]=N
                #triwatersn.append(ii)
                triwatersn[icoast]=ii
                icoast+=1
            #else:
            #    print(" triangulate d1 1",sqrt(d1),sqrt(d2),sqrt(d3),sqrt(ddiagnol))
            ii+=1

        nrange = range(icoast,NODES.shape[0])
        self.triwater = NODES[triwatersn[0:icoast]]
        #self.triwater = np.delete(self.triwater,nrange,0)
        
        if DEBUG :
            print("triwatersn",icoast,triwatersn[icoast])
            print (" triwater     ", self.triwater.shape, self.triwater[:3][1:3])
            print (" tri.simplices", self.tri.simplices.shape)
            print (" self.a", self.a.shape)


# multiprocess:
def buildmesh(meshestobuild, meshesfinished):
    #while True:
        #try:
            # Try to get a task off queue, if queue is empty throw exception and break
            xmm=meshestobuild.get()
        #except queue.Empty:
        #    print(" queue empty",queue.Empty)
        #    break
        #else:
            # do the task using "class DATA(object)  TT" which got popped off the queue

            print("reshape_mesh",xmm.DataType)
            #xmm.get_data(ncl.nc)
            xmm.reshape_mesh()
            xmm.triangulate()
            print("thread mesh done ",xmm.DataType)
            meshesfinished.put(xmm)
            print("thread mesh really done ",xmm.DataType)
            return True


def array_queue(ncfilename,ModelType):
    # set netcdf file name
    # open netcdf file
    # ncfilename='/media/tom/MyBookAllLinux/NOSnetcdf/201902/nos.cbofs.regulargrid.20190201.t00z.n001.nc'
    nc=Dataset(ncfilename)

    if ModelType=="ROMS_REGULAR":
        # Read only, lon, lat, mask
        lon = nc.variables["Longitude"][:]
        lat = nc.variables["Latitude"][:]
        mask= nc.variables["mask"][:]
        xmesh=xyarray(lon,lat,mask,'U',ModelType)
        ymesh=xyarray(lon,lat,mask,'V',ModelType)
        amesh=xyarray(lon,lat,mask,'ANGLE',ModelType)
        wmesh=xyarray(lon,lat,mask,'W',ModelType)
        print(" ROMS_REGULAR lon, lat shapes ",lon.shape,lat.shape)
    elif ModelType=="ROMS_FIELDS":
        lon_u = nc.variables["lon_u"][:]
        lat_u = nc.variables["lat_u"][:]
        mask_u = nc.variables["mask_u"][:]
        xmesh=xyarray(lon_u,lat_u,mask_u,'U',ModelType)
        lon_v = nc.variables["lon_v"][:]
        lat_v = nc.variables["lat_v"][:]
        mask_v = nc.variables["mask_v"][:]
        ymesh=xyarray(lon_v,lat_v,mask_v,'V',ModelType)
        lon_rho = nc.variables["lon_rho"][:]
        lat_rho = nc.variables["lat_rho"][:]
        mask_rho = nc.variables["mask_rho"][:]
        amesh=xyarray(lon_rho,lat_rho,mask_rho,'ANGLE',ModelType)
        wmesh=xyarray(lon_rho,lat_rho,mask_rho,'W',ModelType)
        print(" ROMS_FIELDS lon_u, lon_v shapes ",lon_u.shape,lon_v.shape)

    nc.close()

    # run mask_fix
    # do processes to reshape self.lon self.lat  self.nodes
    # triangulate  add self. tri, a,b,c
    timefirst=time.time()

    # Add the data constructs to the Queue meshestobuild
    meshestobuild=Queue()
    meshesfinished=Queue()
    processes = []

    meshestobuild.put(xmesh)
    meshestobuild.put(ymesh)
    meshestobuild.put(amesh)
    meshestobuild.put(wmesh)

    print(" queue empty?",queue.Empty)

    number_of_processes = 4

    # creating processes
    # Step through the queue pulling out number_of_processes instances:
    # Shouldn't have to count them, but I couldn't make that work
    # So precision counting of processes and puts pull from Queue necessary
    for w in range(number_of_processes):
        p = Process(target=buildmesh, args=(meshestobuild, meshesfinished))
        processes.append(p)
        p.start()
        print("process started ",w)

    for p in processes:
        #p.start()
        print("p.start loop")

    # completing process
    for p in processes:
        print("p.join1")
        #p.join()
        print("p.join2")
    
    #time.sleep(1)    
    print(" done with joins, final state:")
    # print the output

    # Pull the finished data off the Queue, allocate it to 
    #while not meshesfinished.empty():
    for w in range(number_of_processes):
        MM=meshesfinished.get()
        print(MM.DataType,MM.dx)
        if   MM.DataType == 'U':
            xmesh = MM
        elif MM.DataType == 'V':
            ymesh = MM
        elif MM.DataType == 'ANGLE':
            amesh = MM
        elif MM.DataType == 'W':
            wmesh = MM
        else:
            M3 = MM
    
    print("end",time.time()-timefirst)
    return xmesh,ymesh,amesh,wmesh

# Conclusion:
#   The p.join only does harm.
#   meshesfinsihed.get() seems to wait for something to appear
#  If I manually count the number_of_processes then it can work smoothly.

# if this array_queue.py is run by itself this IF executes.
#  If it is included in another mainline this IF should not execute
if __name__ == '__main__':

    ncfilename='/media/tom/MyBookAllLinux/NOSnetcdf/201902/nos.cbofs.regulargrid.20190201.t00z.n001.nc'
    ncfilename='/media/tom/MyBookAllLinux/NOSnetcdf/201905/nos.cbofs.fields.20190513.t18z.n005.nc'
    
    xmesh,ymesh,amesh,wmesh=array_queue(ncfilename,"ROMS_FIELDS")
    
    import matplotlib.pyplot as plt
    import mesh_graphics3 as MG      # Bunch of commented out test plots
    MG.plot_initialize()
    MG.plot_mesh(xmesh,ymesh)
    
    #p,pcolors,numparticles = MG.plot_mesh_pick_line2(xmesh,100)


