# part_mover.py
# routines to move the particles using U24, V24, xmesh, ymesh, amesh, p
# Designed to use multi processing

from multiprocessing import Lock, Process, Queue, current_process
import queue
import multiprocessing as mp

from math import cos, sin, atan2, sqrt
import datetime
import numpy as np


# Split p up over 8 processors
# Send p8[i] to update_part_position
# loop on time
#     itime0, itime1, timefrac
# loop on p8[i][:]
#   find triangle for particle pp
#   interpolate in time the three corners
#   Uinterpolate[0:2] = U[itime0]*(1.-frac) + U[itime1]*frac
#   apply abc's  to get U(p)
#   update p

def loop_updates(xmesh,ymesh,amesh,U24,V24,Time24,p,dt,numtimescalc,numtimesprint,dpinterval):

    timenow=float(Time24[0])
    print("timenow ={0}, Time24={1}  dt={2}".format(timenow,Time24[0],dt))
    itime0=0
    itime1=1
    timefrac=0.0
    timep=[]
    iprint=0
    pp=np.zeros((numtimesprint,np.shape(p)[0],2))
    for idt in range(numtimescalc):

        update_particles(dt,p,xmesh,ymesh,amesh,U24,V24,itime0,itime1,timefrac)
        #print(" idt, p[0] ",p[1])
        
        if idt%dpinterval==0:
            timep.append(timenow)
            pp[iprint]=p
            iprint+=1
        
        # update times and arrays U0, U1,....
        timenow=timenow+dt/(24.*60.*60.)  # units of days, dt in seconds
        timefrac=(timenow-Time24[itime0])/(Time24[itime1]-Time24[itime0])
        if timefrac>1.0 :
            if itime1==len(Time24): break
            itime0+=1
            itime1+=1
            timenow=Time24[itime0]
            timefrac=(timenow-Time24[itime0])/(Time24[itime1]-Time24[itime0])

    return pp,timep
            

def update_particles(dt,p,xmesh,ymesh,amesh,U24,V24,itime0,itime1,timefrac):

    ptrix = xmesh.tri.find_simplex(p)  # find all tri at once
    ptriy = ymesh.tri.find_simplex(p)  # find all tri at once
    ptria = amesh.tri.find_simplex(p)  # find all tri at once
    #print('shape p',np.shape(p))
    #print('ptrix',np.shape(ptrix),ptrix)
    #print('ptriy',np.shape(ptriy),ptriy)
    #print('ptria',np.shape(ptria),ptria)

    ii=0
    for pt in ptrix:
        ptx=ptrix[ii]
        pty=ptriy[ii]
        pta=ptria[ii]
        if ptx>len(xmesh.a):
            ptx=-1
        if pty>len(ymesh.a):
            #print("ii,pty  ",ii,pty,"\n ptriy",ptriy)
            pty=-1
        if pta>len(amesh.a):
            #print("ii,pta  ",ii,pta,"\n ptriy",ptria)
            ptx=-1
            
        #print("pt ",pt)
        #pt=ptri[ii]  # triangle p[ii,0,1] is inside of
        if (ptx>=0 and pty>=0) :    # inside grid
            f0=xmesh.a[pt,0]*p[ii,0] + xmesh.b[pt,0]*p[ii,1] +xmesh.c[pt,0]
            f1=xmesh.a[pt,1]*p[ii,0] + xmesh.b[pt,1]*p[ii,1] +xmesh.c[pt,1]
            f2=xmesh.a[pt,2]*p[ii,0] + xmesh.b[pt,2]*p[ii,1] +xmesh.c[pt,2]
            nodetri=xmesh.tri.simplices[pt,:]     #  three node indices for the triangle
            U0=f0*U24[itime0,nodetri[0]] + f1*U24[itime0,nodetri[1]]+f2*U24[itime0,nodetri[2]]
            U1=f0*U24[itime1,nodetri[0]] + f1*U24[itime1,nodetri[1]]+f2*U24[itime1,nodetri[2]]
            Up=U0*(1.-timefrac) + U1*timefrac
            if ymesh.ModelType=='ROMS_REGULAR':
                V0=f0*V24[itime0,nodetri[0]] + f1*V24[itime0,nodetri[1]]+f2*V24[itime0,nodetri[2]]
                V1=f0*V24[itime1,nodetri[0]] + f1*V24[itime1,nodetri[1]]+f2*V24[itime1,nodetri[2]]
                Vp=V0*(1.-timefrac) + V1*timefrac
#  No rotation for ROMS_REGULAR
                a0=0.0
                Ur=Up
                Vr=Vp
   
            elif ymesh.ModelType=='ROMS_FIELDS':
                f0=ymesh.a[pty,0]*p[ii,0] + ymesh.b[pty,0]*p[ii,1] +ymesh.c[pty,0]
                f1=ymesh.a[pty,1]*p[ii,0] + ymesh.b[pty,1]*p[ii,1] +ymesh.c[pty,1]
                f2=ymesh.a[pty,2]*p[ii,0] + ymesh.b[pty,2]*p[ii,1] +ymesh.c[pty,2]
                nodetri=ymesh.tri.simplices[pty,:]     #  three node indices for the triangle
                V0=f0*V24[itime0,nodetri[0]] + f1*V24[itime0,nodetri[1]]+f2*V24[itime0,nodetri[2]]
                V1=f0*V24[itime1,nodetri[0]] + f1*V24[itime1,nodetri[1]]+f2*V24[itime1,nodetri[2]]
                Vp=V0*(1.-timefrac) + V1*timefrac
                #Urot, Vrot = amesh.rotateangle(p,U,V)   # when needed....
                f0=amesh.a[pta,0]*p[ii,0] + amesh.b[pta,0]*p[ii,1] +amesh.c[pta,0]
                f1=amesh.a[pta,1]*p[ii,0] + amesh.b[pta,1]*p[ii,1] +amesh.c[pta,1]
                f2=amesh.a[pta,2]*p[ii,0] + amesh.b[pta,2]*p[ii,1] +amesh.c[pta,2]
                nodetri=amesh.tri.simplices[pta,:]     #  three node indices for the triangle
                a0=f0*amesh.angle[nodetri[0]] + f1*amesh.angle[nodetri[1]]+f2*amesh.angle[nodetri[2]]

                # a0 is angle between XI axis and East, radians
                ca=cos(a0+4.*3.14159/2.)
                sa=sin(a0+4.*3.14159/2.)
                Ur=ca*Up +sa*Vp
                Vr=-sa*Up +ca*Vp

#            meterperlat = 111132.92 -559.82*cos(2.*p[ii,1].mean())+1.175*cos(4.*p[ii,1].mean())
#            meterperlat = 111132.92 -559.82*cos(2.*p[ii,1].mean())+1.175*cos(4.*p[ii,1].mean())
            meterperlon = 111132.92
            meterperlat = meterperlon*cos(3.14159*p[ii,1]/180.)            
            p[ii,0] = p[ii,0] + dt*Ur/meterperlat
            p[ii,1] = p[ii,1] + dt*Vr/meterperlon
            # Beached particle test, move them out off the Earth
            # Beached particle will have ptx=-1
            if (abs(Ur)+abs(Vr))<.00001 :
                print("move out of space",ii,p[ii])
                p[ii]=(-75.9,37.2)
        else:
            Up=0.0
            Vp=0.0    # outside grid, no triangle, beached
            p[ii]=(-75.8+float(ii)/300.,36.05)
            #p[ii]=(-76.3+float(ii)/300.,37.2)

        ii+=1
    return


# Parallel processing

class pp_object(object):
    
    def __init__(self,p,pcolors,npf,npl,xmesh,ymesh,amesh,U24,V24,Time24,dt,numtimescalc,numtimesprint,dpinterval):
        # These two arrays will be reassembled when processes join up
        self.pp=[]     # np.ones((numtimesprint,np.shape(p)))
        self.timep=[]
        self.pcolors=pcolors
        print("pcolors {0},  self.pcolors{1}".format(np.shape(pcolors),np.shape(self.pcolors)))
        # all these others are just to pass data to the update_particles routine
        # or not
        self.p = p
        self.xmesh=xmesh
        self.ymesh=ymesh
        self.amesh=amesh
        self.U24=U24
        self.V24=V24
        self.Time24=Time24
        self.dt=dt
        self.numtimescalc=numtimescalc
        self.numtimesprint=numtimesprint
        self.dpinterval=dpinterval
        self.npf=npf
        self.npl=npl
                        
def step_particles_forward(pptobuild,ppfinished):

    pobj=pptobuild.get()
    
    pobj.pp,pobj.timep=loop_updates(pobj.xmesh,pobj.ymesh,pobj.amesh,pobj.U24,pobj.V24,pobj.Time24,
                                    pobj.p,pobj.dt,pobj.numtimescalc,pobj.numtimesprint,pobj.dpinterval)

    ppfinished.put(pobj)

    print("thread mesh done ")
    

                        
def queue_pp_calculationloop(p,pcolors,xmesh,ymesh,amesh,U24,V24,Time24,dt,numtimescalc,numtimesprint,dpinterval ):
    # Split p into num_processors and pcolors
    # Create a Queue for pp
    # call loop_updates(xmesh,U24,V24,Time24,p,dt,numtimescalc,numtimesprint,dpinterval)
    # which will put pp to queue

    pptobuild=Queue()
    ppfinished=Queue()
    processes=[]


    # split p in pp_object 's
    number_of_processes = 9
    number_of_particles=np.shape(p)[0]
    npf=0
    nplstep=int(number_of_particles/(number_of_processes-1))
    npl=nplstep
    partspp=[]
    
    for nproc in range(number_of_processes):
        #pp=p[(numtimesprint,range(npf,npl))]   # pp[ntime,nparts,x,y]
        pc = pcolors[range(npf,npl)]
        ppart = p[range(npf,npl)]
        npf=npl
        npl=int(min(npl+nplstep,(number_of_particles)))
        print("nproc={0}, npf={1}, npl={2}, shape(p)={3}, ntime={4}".format(nproc,npf,npl,np.shape(p), numtimesprint))

        pobj = pp_object(ppart,pc,npf,npl,xmesh,ymesh,amesh,U24,V24,Time24,dt,numtimescalc,numtimesprint,dpinterval)
        pptobuild.put(pobj)

        partspp.append(pobj)

        proc_pp=Process(target=step_particles_forward,args=(pptobuild,ppfinished))
        processes.append(proc_pp)
        proc_pp.start()
        
    # completing processes  didn't work last time
    #for p in processes:
    #    print("p.join1")
    #    #p.join()
    #    print("p.join2")

    #time.sleep(1)    
    print(" done with process, rejoin Queue:")
    ppp=np.zeros((numtimesprint,number_of_particles,2))
    ppcolors = np.zeros(number_of_particles,"int")
    for w in range(number_of_processes):
        MM=ppfinished.get()
        #print(" Processor={0}, pp=xx, shape(pp)={2}, {3}".format(w,MM.pp[1],np.shape(MM.pp),np.shape(MM.pcolors)))
        iip=0
        for ip in range(MM.npf,MM.npl):
            ppp[:,ip,:]=MM.pp[:,iip,:]
            ppcolors[ip]=MM.pcolors[iip]
            iip+=1
            
    print(" ppp, shape(ppp)={0}, ppcolors{1}, MM.pcolors{2}".format(np.shape(ppp),np.shape(ppcolors),np.shape(MM.pcolors)))

    return ppp,ppcolors,MM.timep
