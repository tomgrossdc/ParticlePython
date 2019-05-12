# Could build a script to do this:
# wget https://opendap.co-ops.nos.noaa.gov/netcdf/cbofs/201807/nos.cbofs.regulargrid.n003.20180722.t18z.nc

# the Web site to visit is https://opendap.co-ops.nos.noaa.gov/index.php?dir=/netcdf/cbofs/201807/
# subprocess http://www.pythonforbeginners.com/os/subprocess-for-system-administrators

# mkdir 201807
# cd 201807

import subprocess
import shutil, os
import time
from pathlib import Path
from subprocess import call
#call(["ls","-1"])

wwwname="https://opendap.co-ops.nos.noaa.gov/netcdf/cbofs/"
filename="nos.cbofs.regulargrid."

for i in [5,4,3,2,1]:
	ti = time.gmtime(time.time()-i*24*60*60)
	dateget="%d%02d%02d."%(ti.tm_year ,ti.tm_mon,ti.tm_mday)
	yearn="%4d%02d/"%(ti.tm_year,ti.tm_mon)
	print ("Time yearn : ", yearn)
	print ("Time dateget : ", dateget)
	print
	yearnd="%4d%02d"%(ti.tm_year,ti.tm_mon)
	dirname="/media/tom/MyBookAllLinux/NOSnetcdf/"+yearn

	#os.makedirs(dirname,exist_ok=True)
	try:
		call(["mkdir",dirname])
	except :
		print("directory already exists")
			
	
	for sixget in ["t00z.", "t06z.", "t12z.", "t18z.", ] :
		for nh in ["n001.","n002.","n003.","n004.","n005.","n006."] :
			#print(type(filen))
			filen=filename+nh+dateget+sixget+"nc"
			filenew=filename+dateget+sixget+nh+"nc"
			testfile=Path(dirname+filenew)
			if testfile.exists() == False:
				#file does not exist, make a new one
				#print("filen="+filen)
				#print ("wwwname="+wwwname)

				calllist=["wget",wwwname+yearn+filen]
				print (calllist)
				call(calllist)

				calllistmv=["mv",filen,dirname+filenew]
				print ( calllistmv)
				call(calllistmv)
			else:
				print("EXISTS ",dirname+filenew," EXISTS")

ti = time.gmtime()
print ("Time ending  is : "+time.ctime())





