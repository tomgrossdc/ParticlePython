# Read Velocity files 
# file_reader.py

import datetime
import numpy as np
import os, sys
from netCDF4 import Dataset

def ReadOptions():
    filename="/home/tom/code/ThreadParticle/Git20190512/ParticlePython/Particle_Options.txt"
    file1 = open(filename,"r")
    # all data is returned as text so many conversions needed
    file1.readline()  # DEBUG
    DEBUG= bool(file1.readline()[:-1]=='True')
    file1.readline()  # 
    s_rho = int(file1.readline()[:-1])
    file1.readline()  # 
    ModelType = file1.readline()[:-1]
    file1.readline()  # 
    date_start=[int(i) for i in file1.readline()[:-1].split(',')]
    file1.readline()  # 
    date_end=[int(i) for i in file1.readline()[:-1].split(',')]
    file1.readline()  # 
    dt = float(file1.readline()[:-1])
    file1.readline()  # 
    dpsave = float(file1.readline()[:-1])
    file1.readline()  # 
    dirroot = file1.readline()[:-1]
    file1.readline()  # 
    ByHand= bool(file1.readline()[:-1]=='True')
    file1.readline()  # 
    lonlatbox=[float(i) for i in file1.readline()[:-1].split(',')]
    file1.readline()  # 
    lonlatgraphicbox=[float(i) for i in file1.readline()[:-1].split(',')]
    file1.readline()  # 
    nxny=[int(i) for i in file1.readline()[:-1].split(',')]
    file1.readline()  # 
    GRAPHrecord= bool(file1.readline()[:-1]=='True')
    
    file1.close()
    return DEBUG,s_rho,ModelType,date_start,date_end,dt,dpsave,dirroot,ByHand,lonlatbox,lonlatgraphicbox,nxny,GRAPHrecord


def BuildFileList(dirroot,ModelType,datestart,dateend):
    # datestart and dateend are objects from datetime
    # datestart=datetime.datetime(2018,11,10,00,1)
    filenames=[ ]
    filedirectories=[]
    if ModelType=="ROMS_REGULAR" :
    #    filedirectory=("/home/tomg/code/ParticleTracker%d/CRONFILES/%d%d"
    #                   %(datestart.year,datestart.year,datestart.month))
        filedirectory=("%s/%d%02d"
                       %(dirroot,datestart.year,datestart.month))
        namesroot="nos.cbofs.regulargrid."

    deltime=dateend-datestart
    print(deltime.days, deltime.seconds)
    numfiles = int((deltime.days*24*3600+deltime.seconds)/3600 +.5 )
    print("numfiles=%d"%(numfiles))
    del6=datetime.timedelta(hours=1)
    dn=datestart
    i=0
    for hourcount in range(0,numfiles):
        filedirectory, filename=BuildFileName(dirroot,namesroot,dn)
        #print(i,filedirectory," ",filename)
        filenames.append(filename)
        filedirectories.append(filedirectory)
        i+=1
        dn=dn+del6

    return filedirectories, filenames

def BuildFileName(dirroot,namesroot,dn):
    midhour = 6*int(dn.hour/6)
    sixhour=dn.hour-midhour+1
    #print(dn,midhour,sixhour)
    filename=("%s%d%02d%02d.t%02dz.n%03d.nc"
    %(namesroot,dn.year,dn.month,dn.day,midhour,sixhour))
    #filedirectory=("/home/tom/code/ParticleTracker%d/CRONFILES/%d%d"
    #                   %(dn.year,dn.year,dn.month))
    filedirectory=("%s/%d%02d"%(dirroot,dn.year,dn.month))
    return filedirectory, filename

def UV24(filedirectories,filenames,fileindex,mesh,s_rho=0,timedata=0):
    # U24,V24 = FR.UV24(filedirectories,filenames,fileindex,meshx)
    # open timeindex filenames and append to U,V
    # ROMS_REGULAR
    # uses only one mesh.mask  (could be meshx, meshy in future)
    # timedata could be used to access time dim of ROMSnetcdf, ROMS_REGULAR: 0
    # s_rho selects vertical level of the data, ROMS_REGULAR: 0
    print ("start UV24")
    lenti=len(fileindex)
    #lenarray=int(np.sum(mesh.mask)+.001)
    lenarray=len(mesh.mask[mesh.mask>0])
    print('lenarray',lenarray)
    Coast=np.argwhere(mesh.mask>5)
    print("Coast",Coast.shape)

    U24=np.zeros((lenti,lenarray))
    V24=np.zeros((lenti,lenarray))
    Time24=np.zeros(lenti)
    it=0
    for itime in fileindex:
        nclfilename=os.path.join(filedirectories[itime],filenames[itime])
        #print("opening:",nclfilename)
        ncl=Dataset(nclfilename)
        speed = ncl.variables["u_eastward"][timedata,s_rho,:,:]
        #print("speed shape",speed.shape,"mask shape",mesh.mask.shape)
        speed = speed.reshape(speed.shape[0]*speed.shape[1])
        speed[Coast]=0.0
        U24[it] = speed[mesh.mask>0]
        speed = ncl.variables["v_northward"][timedata,s_rho,:,:]
        speed = speed.reshape(speed.shape[0]*speed.shape[1])
        speed[Coast]=0.0
        V24[it] = speed[mesh.mask>0]

        # seconds since 2016-01-01 00:00:00
        timenow=ncl.variables["ocean_time"][timedata]
        # days since 2016-01-01 00:00:00
        Time24[it] = timenow/(24.*3600.)
        ncl.close()
        it+=1

    return U24,V24,Time24






    return U,V,dayhours
