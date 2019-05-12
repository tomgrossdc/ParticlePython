# ParticlePython
Netcdf Python Particle tracker with Multiprocessor and Animations
Runs with AMD Ryzen 7 2700x using all 8 processors.
Does particle tracking in Chesapeake Bay using NOAA netcdf files. 

Missing from the project are those data files. They may be downloaded with the 
crontab tab job, whose files are included in the CRONTAB directory.
The ROMS_REGULAR are ROMS runs which have been meshed to fixed box mesh for 
simplicity and are not good for actual modeling work. But are only hourly files I can find.
There are six hourly real ROMS files. Will make sure the programme will work with them for when I
do find some proper hourly files.

python3 main_mesh.py   will run the main program. 
Options are hardwired inside the code. They should be moved to an input file and read at start of main_mesh.py
Need to specify start time, end time, directory with files. 
Either select a small batch of particles ByHand or generates a box of particles in lonlatbox
Specify the lonlatbox and numx numy of particles. 
200, 100 have worked nicely and runs in about two minutes for five days of simulation.
Output is an animation of particles moving, either on screen or the GRAPHrecord=True generates ches.mp4

Options:  
s_rho = 0    Top layer of particles
ModelType="ROMS_REGULAR"
datestart=(2019,2,1,1,0)
dateend= (2019,2,7,12,0)
dt=120.        Time step, seconds
dpsave=1800.   Time step to record
dirroot = "/media/tom/MyBookAllLinux/NOSnetcdf"
ByHand = True
  Interactive screen position of particles
  ==False
  lonlatbox=(-76.4,36.8,-75.9,37.5)    # Big Box near Mouth
  nx, ny = 100,200  or maybe 10,20 for faster runs
GRAPHrecord = False   plot graphic animation to screen
==True record graphic animation to mp4 file.
