#!/bin/sh
# CronJob.sh
# add to your crontab -e 
# 13 20 * * * /media/tom/MyBookAllLinux/NOSnetcdf/CRONTAB/CronJob.sh &> /media/tom/MyBookAllLinux/NOSnetcdf/CRONTAB/cronjoblog.log

export WGETWORKDIR=/media/tom/MyBookAllLinux/NOSnetcdf/CRONTAB
echo " CRONJob.sh with wget_regular_test.py WGETWORKDIR = ", $WGETWORKDIR
cd /media/tom/MyBookAllLinux/NOSnetcdf/CRONTAB
python3 wget_regular_tomMy.py >> /media/tom/MyBookAllLinux/NOSnetcdf/CRONTAB/junk.log
date
echo "THE END"

