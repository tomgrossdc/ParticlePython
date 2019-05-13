# main_mesh.py 

import array_queue as AQ
import mesh_graphics3 as MG      # Bunch of commented out test plots
import file_reader as FR
import particle_mover_multi as PM
#    import mesh_animate_record as MAR   # These are incompatible
#    import mesh_animate3 as MA          # So are loaded by if else below

import datetime, os, sys, time
import numpy as np
import matplotlib.pyplot as plt

#from mesh_graphics3 import plot_velocities, plot_mesh, plot_initialize


L = FR.ReadOptions()
print("L",L)
DEBUG=L[0]
s_rho=L[1]
ModelType=L[2]
date_start =L[3]
date_end =L[4]
dt =L[5]
dpsave =L[6]
dirroot =L[7]
ByHand =L[8]
lonlatbox =L[9]
lonlatgraphicbox =L[10]
nxny =L[11]
GRAPHrecord=L[12]

# from starttime, lasttime  create list of file names

# Layer to pick from the netcdf files
# ROMS_REGULAR is 0-top   15-deepest,
# If depth is less than depth(s_rho), U=0
#s_rho=0

#ModelType="ROMS_REGULAR"
datestart=datetime.datetime(*date_start)
dateend = datetime.datetime(*date_end)
timespan=dateend-datestart             # datetime construct for differences
timespansec=timespan.total_seconds()   #  real seconds
#dt = 120.  # seconds per calculation time step
#dpsave = 1800. # seconds per particle position save
numtimescalc=int(timespansec/dt)
numtimesprint=int(timespansec/dpsave)
dpinterval = int(dpsave/dt +.01)
print(" Timings:\n Duration(sec)={0}\n dt={1}  dpsave={2}\n numtimescalc={3}, numtimesprint={4}, dpinterval={5}"
      .format(timespansec,dt,dpsave,numtimescalc,numtimesprint,dpinterval))


      
#dirroot="/media/tom/MyBookAllLinux/NOSnetcdf"
filedirectories, filenames=FR.BuildFileList(dirroot,ModelType,datestart,dateend)
print("BuildFileList says",filedirectories[0])
print("BuildFileList says",filenames[0])

nclfilename=os.path.join(filedirectories[0],filenames[0])

# Read the first netcdf file, pull out mesh data, triangulate and prepare meshes
xmesh,ymesh,amesh=AQ.array_queue(nclfilename)

MG.plot_initialize()

#p = MG.plot_mesh_pick_line(xmesh,numtimes)
#ByHand=False
if ByHand:
    p,numparticles = MG.plot_mesh_pick_line2(xmesh,numtimesprint)
    print("p {0}, numparticles={1}\n".format(np.shape(p),numparticles),p[:,:])
else:
    #lonlatbox=(-76.04,36.95,-76.005,37.15)  # Just inside Mouth slender tall
    #lonlatbox=(-76.20,36.95,-76.005,37.15)  # Just inside Mouth big box
    #lonlatbox=(-76.16,37.59,-76.08,37.64)
    #lonlatbox=(-76.46,38.67,-76.40,38.73)
    #lonlatbox=(-76.51,38.67,-76.36,38.73)   # long box near annapolis
    #lonlatbox=(-76.43,38.70,-76.40,38.73)
    #lonlatbox=(-76.6,38.2,-76.1,39.2)    # Too Big box near Potomac
    #lonlatbox=(-76.4,36.8,-75.9,37.5)    # Too Big Box near Mouth
    #lonlatbox=(-76.55,36.8,-75.9,39.2)    # Really Big Box from Mouth to Potomac
    p,pcolors,numparticles=MG.plot_box_pick(nxny[0],nxny[1],lonlatbox) 
    #p,pcolors,numparticles=MG.plot_box_pick(100,200,lonlatbox) # 164 sec 2000x216X1800sec


# initialize first U,V by openning first two files
#open netcdf for U,V
# U1,time0 = xmesh.velocity(ncl.nc)
# U = xmesh.findvel(p,U1)
# Urot, Vrot = amesh.rotateangle(p,U,V)


# Parallelize on time or particles
# open 24 hours of U,V,using netcdf and put into big constructs
# ALLMESHES= UV24(filedirs,filenames,timeindex)
#launch batches of particles
#  simple np.array Particles[NUMTime][NUMPart][X, .Y, .Z ]
# Direct to save positions only every hour, but dt=120sec
# fileindex=range(0,25)  #  gives 0:24;  25 numbers
# fileindex=range(24,49) # gives 24:48; read 25 to give overlap with previous group
# This is fast and four days barely dents memory, so grab all data first.
# In future add method to read more during a run.


timefirst=time.time()
ifilemax=len(filenames) 
fileindex=range(0,ifilemax)
U24,V24,Time24 = FR.UV24(filedirectories,filenames,fileindex,xmesh,s_rho)

"""import matplotlib.pyplot as plt
plt.plot(U24[4])
plt.show
"""

print("time for max",ifilemax,'time',time.time()-timefirst)
print('U24',np.shape(U24))
print('Time24',Time24)
print('U24,V24[3][250]',U24[3][250],V24[3][250])

# lousey debug plot
#MG.plot_velocities(xmesh,U24[4],V24[4],p,U24[4],V24[4])

# Now have U24,V24[:itime,:NumNodes]
# p[:numparticles, :0,1 X,Y]
print(' Loop {0} times \n particles={1}'.format(numtimescalc,np.shape(p)))


t1=time.time()
pp,ppcolors,timep= PM.queue_pp_calculationloop(p,pcolors,xmesh,U24,V24,Time24,dt,numtimescalc,numtimesprint,dpinterval)
print('Parallel time loop updates',time.time()-t1, "seconds")

#print('shape(pp)={0},\n pp={1}\n'.format(np.shape(pp),pp))


#MG.plot_particles(xmesh,pp,pcolors)

#MG.plot_particles_spots(xmesh,pp,pcolors)
#print("timep[3]-",timep[3])

# lonlatgraphicbox=(-76.6,36.2,-75.6,37.6)
if GRAPHrecord :
    import mesh_animate_record as MAR
    manim=MAR.mesh_animate(lonlatgraphicbox)  # (-76.6,36.2,-75.6,37.6)
    #manim.get_particles_time()
    if DEBUG:
        snps = len(pcolors)
        sntime = int(len(pp)/snps)
        print("get_particles_pp, n_particles={0}, n_time={1}, pp={2}".format(snps, sntime,np.shape(pp)))
    manim.get_particles_pp(pp,pcolors,timep)
    manim.run_animation(xmesh,100)
    #manim.run_animation(False,100)
else:
    import mesh_animate3 as MA
    manim=MA.mesh_animate(lonlatgraphicbox)   # (-76.6,36.4,-75.5,38.0)
    #manim.get_particles_time()
    if DEBUG:
        snps = len(pcolors)
        sntime = int(len(pp)/snps)
        print("get_particles_pp, n_particles={0}, n_time={1}, pp={2}".format(snps, sntime,np.shape(pp)))
    manim.get_particles_pp(pp,pcolors,timep)
    manim.run_animation(xmesh,100)
    #manim.run_animation(False,100)
