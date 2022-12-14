import os,sys
import glob
import numpy as np
import matplotlib,netCDF4
import matplotlib.pyplot as plt
import datetime
import struct
import pandas as pd
import math
import scipy
from scipy.sparse.linalg import lsqr
from scipy.spatial.transform import Rotation as R
import time as timeit
import shutil
from rtp_adcp_class import *

## Setting time origin
rtime=datetime.datetime(2020,1,1,0,0,0)



## Set path to save data - netcdf files to give back to glider
odir='./files_FOR_glider/'

# Set path to PDO files from glider
idir = './files_FROM_glider/'

## Set path for completed files
pdir = './files_processed/'

#grab most recent file
files=glob.glob(idir+'*.pd0')
files.sort(key=os.path.getmtime)


all_run_time = []
for f in files: 
    run_time_s = timeit.time()
    file = f 
    
    if (os.path.isfile(file)) & (os.stat(file).st_size>0):
    

        #to access variables outside of class, must assign variable to initialization also need it to access functions like rdi_proc.readPD0
        #in this case rdi_proc, so if I want to access a variable outside of the class (i.e. anything inside of self)

        rdi_proc=RDI_real_time_proc(file,idir,odir)
        rdi_proc.read_PD0(file)


        # set variables for inversion
        U=rdi_proc.u1.transpose() #u from read_PD0 output
        V=rdi_proc.u2.transpose() #v from read_PD0 output
        dz=10
        u_daverage=0
        v_daverage=0
        bins=rdi_proc.bins #bins from read_PD0 output
        depth=rdi_proc.depth #depth from read_PD0 output
        wDAC=5
        wSmoothness=1


        # check that depth is not stuck at the surface and check is data is all nans and percentage of nans
        check_depth = rdi_proc.depth
        check_u = U

        err_msg,err_flag = misc_checks(check_depth, check_u)


        if err_flag == 1 or err_flag ==2:
            print(err_msg)
        else:
            print(err_msg)



            #initialize class
            glider_inv = Glider_ADCP_Inversion(U,V,dz,u_daverage,v_daverage,bins,depth, wDAC, wSmoothness)
            #run inversion
            glider_inv.inversion()

            #define vars for writing netcdf
            O_ls = glider_inv.O_ls
            G_ls = glider_inv.G_ls
            bin_new = glider_inv.bin_new
            time = rdi_proc.time
            #initialize class
            write_nc = Write_NC(file,odir,O_ls, G_ls, bin_new,time)
            #create netcdf
            write_nc.write_data()

            #define vars for plotting
#             depth = rdi_proc.depth
#             pitch = rdi_proc.pitch
#             roll = rdi_proc.roll
#             heading = rdi_proc.heading
#             bins = rdi_proc.bins
#             u1=rdi_proc.u1
#             u2=rdi_proc.u2
#             u3=rdi_proc.u3
#             u4=rdi_proc.u4
#             c1=rdi_proc.c1
#             c2=rdi_proc.c2
#             c3=rdi_proc.c3
#             c4=rdi_proc.c4
#             ei1=rdi_proc.ei1
#             ei2=rdi_proc.ei2
#             ei3=rdi_proc.ei3
#             ei4=rdi_proc.ei4

#             plot_dvl = Plot_ADCP(time, pitch, roll, heading, depth, bins,u1,u2,u3,u4,c1,c2,c3,c4,ei1,ei2,ei3,ei4,O_ls,bin_new)
#             plot_dvl.plot_data()
            
            #move processed file into processed folder
            shutil.move(f,pdir)
            
            run_time_e = timeit.time()

            time_tot = run_time_e-run_time_s

            all_run_time=np.append(all_run_time,time_tot)
            print('Finished')
    else:
        print('File does not exist.')


t_list = np.array(all_run_time)



print('min time ',round(t_list.min(),3))

print('max time ',round(t_list.max(),3))

print('avg time ',round(t_list.mean(),3))