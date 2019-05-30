# ParticlePython
Netcdf Python Particle tracker with Multiprocessor and Animations
Runs with AMD Ryzen 7 2700x using all 8 processors.
Does particle tracking in Chesapeake Bay using NOAA netcdf files. 

Missing from the project are those data files. They may be downloaded with the 
crontab tab job, whose files are included in the CRONTAB directory.
The ROMS_REGULAR are ROMS runs which have been meshed to fixed box mesh for 
simplicity and are not good for actual modeling work. 

Recently discovered that ROMS_FIELDS files exist. They are hourly ROMS on the
native ROMS grid.  Good God that grid is bizarre. The ROMS rectangular curvilinear
grids are not basically discontinuous, distorted so that rivers can be bent. 
Triangulation of those grids is a major pain. Solution was to use the mask on the
extra ROMS points, and then add in a bunch of points surrounding the coast. These
are then used to mask the triangles of the Delaunay triangulation. 
Whew. It works now. And it is pretty fast. Five days, 100x600=60,000 particles in 538 seconds.

python3 main_mesh.py   will run the main program. 
Options are read from an input file Particle_Options.txt from file_reader.py
Need to specify start time, end time, directory with files. 
Either select a small batch of particles ByHand or generates a box of particles in lonlatbox
Specify the lonlatbox and numx numy of particles. 
200, 100 have worked nicely and runs in about two minutes for five days of simulation.
Output is an animation of particles moving, either on screen or the GRAPHrecord=True generates ches.mp4

Options:  
#DEBUG 
True
#s_rho 0=top:Regular  19=top:Fields
19
#ModelType  ROMS_REGULAR or ROMS_FIELDS
ROMS_FIELDS
#date_start
2019,5,23,7,0
#date_end
2019,5,25,19,0
#dt
120.
#dpsave
1800.
#dirroot
/media/tom/MyBookAllLinux/NOSnetcdf
#ByHand   
False
#lonlatbox Really Big Box from Mouth to Potomac
-76.55,36.8,-75.9,39.2
#lonlatgraphicbox
-76.8,36.5,-75.8,39.5
#nxny
25,50
#GRAPHrecord
True

