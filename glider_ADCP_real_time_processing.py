## Imports
#!/usr/bin/env python
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

## Setting time origin
rtime=datetime.datetime(2020,1,1,0,0,0)

## Set path to data
odir='./'
idir='./'



time=[]    
depth=[] 
pitch=[]
roll=[]
heading=[]
temp=[]      
u1=[]
u2=[]
u3=[]
u4=[] 
ei1=[]
ei2=[]
ei3=[]
ei4=[] 
c1=[]
c2=[]
c3=[]
c4=[]  
pg1=[]
pg2=[]
pg3=[]
pg4=[]       
bins=[]
O_ls=[]
G_ls=[]          
plt.ion()


## Main function to run through entire processing workflow
def main(argv):
    files=glob.glob(idir+'*.pd0')
    files.sort(key=os.path.getmtime)
    files=files[0]
    
    #for file in files:
    read_PD0(files)
    #process_data(U=u1,V=u2,H=1000,dz=10,u_daverage=0,v_daverage=0)
    inversion(U=u1.transpose(), V=u2.transpose(), dz=10, u_daverage=0, v_daverage=0, bins=bins, depth=depth, wDAC=5, wSmoothness=1)
    write_data(files)  
    plot_data()
    plt.show()


## Read in binary PD0 file
def read_PD0(infile):
    global time,depth,pitch,roll,heading,temp,bins
    global u1,u2,u3,u4
    global c1,c2,c3,c4
    global ei1,ei2,ei3,ei4
    global pg1,pg2,pg3,pg4
    print('Reading PDO file : '+infile)
    
    ## Open file and read in binary
    f=open(infile,'rb')
    dat = f.read()
    f.close()
    
    ## All this does is try to find the byte location of the first ensemble. 
    ## Usually its the first byte, but it some applications it is not.  
    ## It searches for the header ID and source ID (both '0x7f')
    ## See Chapter 8 of  PathFinder DVL Guide_Apr20.pdf 
    for [ind,tmp1] in enumerate(dat):
          if  (hex(dat[ind])=='0x7f') and (hex(dat[ind+1]) =='0x7f'):                 
                  break
    
    ## This extracts the number of bytes per ensemble. 
    nbytes=struct.unpack("h",dat[ind+2:ind+4])[0]+2
#    print('Finding Ensembles')      

    ## Find the starting byte of every ensemble. (Varaible Iens)
    ## It goes through every byte and searchs for header ID's and source ID's
    ## When one is found, the index is added to the variable Iens 
    Iens=[] ## Starting byte of each ensemble. 
    n=0

    for [ind,tmp1] in enumerate(dat):

          if ind == len(dat)-1: break

          if (hex(dat[ind])=='0x7f') and (hex(dat[ind+1]) =='0x7f'):
              n=n+1
              nbytes2=struct.unpack("h",dat[ind+2:ind+4])[0]+2  
              startens=ind
              tdat=dat[startens:startens+nbytes]  
              if len(tdat)<nbytes:
                   print('breaking')
                   break
              tmp=tdat[nbytes-2:nbytes]
              chksum=struct.unpack("<H",tmp)[0]
              if (sum(tdat[:nbytes-2]) & 65535) ==  chksum:
                      
                   if nbytes == nbytes2:
  

                       Iens.append(ind)
 #             else:
 #                print('Bad Checksum')

    nens=len(Iens)

## nens is number of ensembles, so this chunk is preallocating variables that will
## be read in. Is the 100 for the 2-dimensional variables just a "safe" number of bins
## without knowing exactly how many. There's a better way to do this...
    time=np.empty((nens),np.double)    
    depth=np.empty((nens),np.double) 
    pitch=np.empty((nens),np.double) 
    roll=np.empty((nens),np.double) 
    heading=np.empty((nens),np.double) 
    temp=np.empty((nens),np.double)
    
    u1=np.empty((nens,100),np.double) 
    u2=np.empty((nens,100),np.double) 
    u3=np.empty((nens,100),np.double)
    u4=np.empty((nens,100),np.double)  
    ei1=np.empty((nens,100),np.double) 
    ei2=np.empty((nens,100),np.double) 
    ei3=np.empty((nens,100),np.double)
    ei4=np.empty((nens,100),np.double)  
    c1=np.empty((nens,100),np.double) 
    c2=np.empty((nens,100),np.double) 
    c3=np.empty((nens,100),np.double)
    c4=np.empty((nens,100),np.double)   
    pg1=np.empty((nens,100),np.double) 
    pg2=np.empty((nens,100),np.double) 
    pg3=np.empty((nens,100),np.double)
    pg4=np.empty((nens,100),np.double)    
    xform=np.zeros((4,4),np.double)   
    xformR=np.zeros((3,3),np.double)
    xformP=np.zeros((3,3),np.double)
    xformH=np.zeros((3,3),np.double)

    ind=0
    eoffset=0
#    Iens=Iens[0:nens]

## Loop through ensembles and pull out data.
## Which bytes correspiond to what variables is detailed in the Pathfinder manual Ch 8. This is standard for PD0s?
    for ind2 in Iens:
        startens=(ind2)
        tdat=dat[startens:startens+nbytes]
        # a=buffer(tdat,2,2)
        tnbytes=struct.unpack("H",tdat[2:4])[0]+2
        # a=buffer(tdat,nbytes-2,2)
        chksum=struct.unpack("<H",tdat[nbytes-2:nbytes])[0]

        ## In the past binary data was subject to erros (missing bytes)during read/writing/data transfer operations. 
        ## It is common practive to add a Checksum to the end of the data. 
        ## i.e. the last 2 bytes of the ensemble represent the "sum" of every byte in the ensemble. 
        ## if they do not match, data was lost in the ensemble. This happens rarely.
        if (sum(tdat[:nbytes-2]) & 65535) ==  chksum:
              ndtype=struct.unpack("b",tdat[5:6])[0]
 
              offsets=list()
              for ind3 in range(ndtype):
                   Is=6+ind3*2
                   offsets.append(struct.unpack_from("h",tdat[Is:Is+2])[0])        
            
            ## FIXEDLEADER
            ## Number of beams
              Is=offsets[0]+8
              nbeam=tdat[Is]
            ## Number of cells
              Is=offsets[0]+9
              ncells=tdat[Is]  
            ## Cell size
              Is=offsets[0]+12
              cellsize=struct.unpack("H",tdat[Is:Is+2])[0]        
              cellsize=cellsize/100.0       
            ## Bin 1 distance in cm --> meters
              Is=offsets[0]+32
              bin1=struct.unpack("H",tdat[Is:Is+2])[0]    
              bin1=bin1/100.0
            ## Heading alignment/100 to get degrees    
              Is=offsets[0]+26
              hdalign=struct.unpack("H",tdat[Is:Is+2])[0]    
              hdalign=hdalign/100.0
            ## Heading bias/100 to get degrees
              Is=offsets[0]+28
              hdbias=struct.unpack("H",tdat[Is:Is+2])[0]    
              hdbias=hdbias/100.0
            
              Is=offsets[0]+4
              # sysconfig1=bin(tdat[Is])   
              sysconfig1=format(tdat[Is], '#010b')[2:]
              Is=offsets[0]+5
              # sysconfig2=bin(tdat[Is]) 
              sysconfig2=format(tdat[Is], '#010b')[2:]
            
            ## Beam angle correction
              if sysconfig2[-2:]=='10':
                  bmang=30.0
              elif sysconfig2[-2:]=='01':
                  bmang=20.0
              elif sysconfig2[-2:]=='00':
                  bmang=15.0
              a=1.0/(2.0*np.sin(bmang*np.pi/180.0))
              b=1.0/(4.0*np.cos(bmang*np.pi/180.0))
              c=1.0
              d=a/np.sqrt(2.0)
            ## Building transformation matrix for beam to instrument
              xform[0,0]=c*a  
              xform[0,1]=-c*a
              xform[0,2]=0.0
              xform[0,3]=0.0  
            
              xform[1,0]=0.0  
              xform[1,1]=0.0
              xform[1,2]=-c*a
              xform[1,3]=c*a
              
              xform[2,0]=b  
              xform[2,1]=b
              xform[2,2]=b
              xform[2,3]=b
              
              xform[3,0]=d  
              xform[3,1]=d
              xform[3,2]=-d
              xform[3,3]=-d

              Is=offsets[1]+2
              ens=struct.unpack("H",tdat[Is:Is+2])[0]        
              
              
     
            ## Read in time data
              Is=offsets[1]+4
              year=tdat[Is]  
        
              Is=offsets[1]+5
              month=tdat[Is]
              Is=offsets[1]+6
              day=tdat[Is]
              Is=offsets[1]+7
              hour=tdat[Is]
              Is=offsets[1]+8
              minute=tdat[Is]
              Is=offsets[1]+9
              sec=tdat[Is]
              Is=offsets[1]+10
              hsec=tdat[Is]

              ttime = datetime.datetime(year+2000,month,day,hour,minute,sec,hsec*10)-rtime
            ## Read in depth data
              Is=offsets[1]+16
              tdepth=struct.unpack("H",tdat[Is:Is+2])[0]   
            ## Convert depth from decimeters to meters
              tdepth=tdepth*0.1
            ## Read in temporary heading
              Is=offsets[1]+18
              theading=struct.unpack("H",tdat[Is:Is+2])[0]     
              theading=theading/100.0
            ## Read in pitch data
              Is=offsets[1]+20
              tpitch=struct.unpack("h",tdat[Is:Is+2])[0]        
              tpitch=tpitch/100.0 
            ## Read in roll data
              Is=offsets[1]+22
              troll=struct.unpack("h",tdat[Is:Is+2])[0]        
              troll=troll/100.0    
            ## Read in salinity data
              Is=offsets[1]+24
              tsalt=struct.unpack("h",tdat[Is:Is+2])[0]
            ## Read in temperature data
              Is=offsets[1]+26
              ttemp=struct.unpack("h",tdat[Is:Is+2])[0]        
              ttemp=ttemp/100.0       
            ## Read in glider pressure data
              Is=offsets[1]+48
              tpress=struct.unpack("i",tdat[Is:Is+4])[0]        
              tpress=tpress/1000.0       
            ## Read in velocity data
              Is=offsets[2]+2
              fmt = "<%dh" % (ncells*4)
              uvw=struct.unpack(fmt,tdat[Is:Is+ncells*4*2])
              uvw=np.array(uvw,dtype=float)
            ## Read in echo intensity data
              Is=offsets[3]+2
              fmt = "<%dB" % (ncells*4)
              tEI=struct.unpack(fmt,tdat[Is:Is+ncells*4])
              tEI=np.array(tEI)
            ## Read in correlation data
              Is=offsets[4]+2
              fmt = "<%dB" % (ncells*4)
              tC=struct.unpack(fmt,tdat[Is:Is+ncells*4])
              tC=np.array(tC)  
            ## Read in percent good data
              Is=offsets[5]+2
              fmt = "<%dB" % (ncells*4)
              tPG=struct.unpack(fmt,tdat[Is:Is+ncells*4])
              tPG=np.array(tPG)              
            ## Reshape velocity, ei, corr, and pg to be 2D based on cells and beams
              uvw.shape=(ncells,4)
              tEI.shape=(ncells,4)
              tC.shape=(ncells,4)
              tPG.shape=(ncells,4)
              
              ## QAQC Data
              uvw = qaqc_data(uvw, tEI, tC, tPG)              
              
            ## Create bins variable based on number of cells and cell size and reference it to distance from sensor based on bin1
              bins=(np.arange(0,ncells,1,np.double)*cellsize)+bin1    
            ## Read in bottom track data
              LO=len(offsets)
            
            ## When the ADCP writes the PD0, each ensemble (output profile) has a set of data types. i.e. fixed leader,
            ## variable leader, velocity, correltation, etc. Each data type has it's own header with an identifier.
            ## When we read a PD0 file, we read each ensemble one at a time.  (variable tdat in the script).
            ## Early in the ensemble, there is a variable called offsets, giving us the byte offset from the beginning
            ## of the  ensemble to the start of the data types in the file. If there is no bottom track, one of these
            ## offsets will be missing. So of LO>6, then there is bottom track  data to process. 

              if LO>6:
                  Is=offsets[6]
                  tmp1=struct.unpack("c",tdat[Is:Is+1])[0]
                  Is=offsets[6]+1
                  tmp2=struct.unpack("c",tdat[Is:Is+1])[0]
                  if [tmp1+tmp2]==[b'\x00\x06']:
                      Is=offsets[6]+16
                      tr1=struct.unpack("H",tdat[Is:Is+2])[0]       
                      Is=offsets[6]+18
                      tr2=struct.unpack("H",tdat[Is:Is+2])[0]    
                      Is=offsets[6]+20
                      tr3=struct.unpack("H",tdat[Is:Is+2])[0]    
                      Is=offsets[6]+22
                      tr4=struct.unpack("H",tdat[Is:Is+2])[0]     
                ## Convert bottom track range to meters
                      tr1=tr1/100.0     
                      tr2=tr2/100.0     
                      tr3=tr3/100.0     
                      tr4=tr4/100.0
                ## Only set velocity data below detected bottom equal to 0 if the bottom track range is greater than 0
                      if tr1>0:
                          uvw[bins>.85*tr1,0]=float("NAN")
                      if tr2>0:
                          uvw[bins>.85*tr2,1]=float("NAN")
                      if tr3>0:
                          uvw[bins>.85*tr3,2]=float("NAN")
                      if tr4>0:
                          uvw[bins>.85*tr4,3]=float("NAN")
                  
            ## Bin map velocity data based on pitch and roll
              uvw=mapdepthcells(uvw,tpitch,troll)
            ## This is where we actually build the arrays/vecotrs of the PD0 data. 
            ## ind is incremented at the end of the code block (ind=ind+1 below)
            ## It's how we add data to the variables we initialized earlier.  
              time[ind]=ttime.days+ttime.seconds/86400.0
                  
              depth[ind]=tdepth
              pitch[ind]=tpitch
              roll[ind]=troll
              temp[ind]=ttemp
              heading[ind]=theading  
              # P=np.arctan(np.tan(tpitch*np.pi/180.0)*np.cos(troll*np.pi/180.0))

            ## Correct heading
              shead=theading+hdalign
            ## Build beam to ENU transformation matrix
              CH=np.cos(shead*np.pi/180.0)
              SH=np.sin(shead*np.pi/180.0)
              CR=np.cos(troll*np.pi/180.0)
              SR=np.sin(troll*np.pi/180.0)
              CP=np.cos(tpitch*np.pi/180.0)
              SP=np.sin(tpitch*np.pi/180.0)
              # print(CP)
              # CP=np.cos(P)
              # print(CP)
              # SP=np.sin(P)
              
              xformH[0,0]=CH  
              xformH[0,1]=SH
              xformH[0,2]=0.0 
              xformH[1,0]=-SH 
              xformH[1,1]=CH
              xformH[1,2]=0.0 
              xformH[2,0]=0.0  
              xformH[2,1]=0.0
              xformH[2,2]=1.0 
            
              xformR[0,0]=CR  
              xformR[0,1]=0.0
              xformR[0,2]=SR 
              xformR[1,0]=0.0 
              xformR[1,1]=1.0
              xformR[1,2]=0.0 
              xformR[2,0]=-SR  
              xformR[2,1]=0.0
              xformR[2,2]=CR
              
              xformP[0,0]=1.0  
              xformP[0,1]=0.0
              xformP[0,2]=0.0 
              xformP[1,0]=0.0 
              xformP[1,1]=CP
              xformP[1,2]=-SP 
              xformP[2,0]=0.0  
              xformP[2,1]=SP
              xformP[2,2]=CP
              
            ## Convert from beam to ENU velocity
              uvw=uvw @ xform.transpose()
              terr=uvw[:,3]
              tuvw=uvw[:,0:3]
              tuvw=tuvw @ xformR.transpose()
              tuvw=tuvw @ xformP.transpose()
              tuvw=tuvw @ xformH.transpose()
            
              # u1[ind,0:ncells]=uvw[:,0]
              # u2[ind,0:ncells]=uvw[:,1]
              # u3[ind,0:ncells]=uvw[:,2]
              # u4[ind,0:ncells]=uvw[:,3]
              u1[ind,0:ncells]=tuvw[:,0] ## This is now U velocity
              u2[ind,0:ncells]=tuvw[:,1] ## This is now V velocity
              u3[ind,0:ncells]=tuvw[:,2] ## This is now W velocity
              u4[ind,0:ncells]=terr      ## This is now error velocity
              ei1[ind,0:ncells]=tEI[:,0] 
              ei2[ind,0:ncells]=tEI[:,1]
              ei3[ind,0:ncells]=tEI[:,2]
              ei4[ind,0:ncells]=tEI[:,3]
              c1[ind,0:ncells]=tC[:,0]
              c2[ind,0:ncells]=tC[:,1]
              c3[ind,0:ncells]=tC[:,2]
              c4[ind,0:ncells]=tC[:,3]
              pg1[ind,0:ncells]=tPG[:,0]
              pg2[ind,0:ncells]=tPG[:,1]
              pg3[ind,0:ncells]=tPG[:,2]
              pg4[ind,0:ncells]=tPG[:,3]            
              
              ind=ind+1
        else:
          #   print 'BAD CHECKSUM'
            # eoffset=eoffset+1
            
            continue
    
    
    u1=u1[:,0:ncells]
    u2=u2[:,0:ncells]
    u3=u3[:,0:ncells]
    u4=u4[:,0:ncells]
    c1=c1[:,0:ncells]
    c2=c2[:,0:ncells]
    c3=c3[:,0:ncells]
    c4=c4[:,0:ncells]
    ei1=ei1[:,0:ncells]
    ei2=ei2[:,0:ncells]
    ei3=ei3[:,0:ncells]
    ei4=ei4[:,0:ncells]
    pg1=pg1[:,0:ncells]
    pg2=pg2[:,0:ncells]
    pg3=pg3[:,0:ncells]
    pg4=pg4[:,0:ncells]
    
    #return(time,depth,pitch,roll,heading,temp,bins,u1,u2,u3,u4,c1,c2,c3,c4,ei1,ei2,ei3,ei4,pg1,pg2,pg3,pg4)


def mapdepthcells(uvw,tpitch,troll):
    global bins
 #   print('Mapping depth cells')
    brange=bins/np.cos(30*np.pi/180.0)
    tuvw=uvw*np.nan
    az=90*np.pi/180
    elev=-60*np.pi/180
    ## This is getting the actuall bin depths.(relative to glider), so we can transform them to "Level" bin depths 
    XYZ1 = sph2cart(az,elev,brange)
    
   
    az=-90*np.pi/180
    elev=-60*np.pi/180
    XYZ2 = sph2cart(az,elev,brange)
    
    az=0*np.pi/180
    elev=-60*np.pi/180
    XYZ3 = sph2cart(az,elev,brange) 
    
    az=180*np.pi/180
    elev=-60*np.pi/180
    XYZ4= sph2cart(az,elev,brange)
#    trot=np.array([[0.9330  , -0.0670 ,  -0.3536],[-0.0670  , 0.9330 ,  -0.3536],[0.3536 ,   0.3536 ,   0.8660]])


#    print([tpitch,troll])
    rang1=tpitch*np.pi/(180) 
    rang2=troll*np.pi/(180) 
    sc=1/np.sqrt(2.0)
    ax1=sc*rang1*np.array([-1, 1, 0])
    ax2=sc*rang2*np.array([1, 1, 0])
    rmx1 = R.from_rotvec(ax1)
    rmx2 = R.from_rotvec(ax2)
    
    rot1=rmx1.as_matrix()    
    rot2=rmx2.as_matrix()
    rXYZ1=np.concatenate(([XYZ1[0]],[XYZ1[1]],[XYZ1[2]]),axis=0)
    r1=rXYZ1.transpose() @ rot1 @ rot2
    rXYZ2=np.concatenate(([XYZ2[0]],[XYZ2[1]],[XYZ2[2]]),axis=0)
    r2=rXYZ2.transpose() @ rot1  @ rot2
    rXYZ3=np.concatenate(([XYZ3[0]],[XYZ3[1]],[XYZ3[2]]),axis=0)
    r3=rXYZ3.transpose() @ rot1 @ rot2
    rXYZ4=np.concatenate(([XYZ4[0]],[XYZ4[1]],[XYZ4[2]]),axis=0)
    r4=rXYZ4.transpose() @ rot1 @ rot2
    
    # gX=np.arange(-10,10)
    # gY=gX*0.0
    # g=(gX+1j*gY)*np.exp(1j*45.0*np.pi/180)
    # gX=np.real(g)
    # gY=np.imag(g)
    # gZ=0.0*gX
    
    
    # gDX=bins*0.0
    # gDY=bins*0.0
    # gDZ=-bins
    # tmp=np.concatenate(([gX],[gY],[gZ]),axis=0)
    # rg=tmp.transpose() @ rot1 @ rot2
    
    # tmp=np.concatenate(([gDX],[gDY],[gDZ]),axis=0)
    # rDg=tmp.transpose() @ rot1 @ rot2
    
    # plt.figure(1)
    # plt.clf()
    # ax = plt.axes(projection='3d')
    # ax.plot3D(XYZ1[0],XYZ1[1],XYZ1[2],'r')
    # ax.plot3D(r1[:,0],r1[:,1],r1[:,2],'b')

    # ax.plot3D(XYZ2[0],XYZ2[1],XYZ2[2],'r')
    # ax.plot3D(r2[:,0],r2[:,1],r2[:,2],'b')
    # ax.plot3D(XYZ3[0],XYZ3[1],XYZ3[2],'r')
    # ax.plot3D(r3[:,0],r3[:,1],r3[:,2],'b')
    # ax.plot3D(XYZ4[0],XYZ4[1],XYZ4[2],'r')
    # ax.plot3D(r4[:,0],r4[:,1],r4[:,2],'b')
    # ax.plot3D(gX,gY,gZ,'k')
    # ax.plot3D(rg[:,0],rg[:,1],rg[:,2],'b')
    

    # ax.plot3D(gDX,gDY,gDZ,'r')
    # ax.plot3D(rDg[:,0],rDg[:,1],rDg[:,2],'b')
    # plt.grid(True)
    
    tuvw[:,0]=np.interp(-r1[:,2],bins,uvw[:,0])
    tuvw[:,1]=np.interp(-r2[:,2],bins,uvw[:,1])
    tuvw[:,2]=np.interp(-r3[:,2],bins,uvw[:,2])
    tuvw[:,3]=np.interp(-r4[:,2],bins,uvw[:,3])
    # plt.figure(3)
    # plt.clf()
    # plt.subplot(141)
    # plt.plot(uvw[:,0],-bins,'r',marker='o')
    # plt.plot(tuvw[:,0],-bins,'r',marker='+')
    # plt.grid(True)
    # plt.subplot(142)
    # plt.plot(uvw[:,1],-bins,'b',marker='o')
    # plt.plot(tuvw[:,1],-bins,'b',marker='+')
    # plt.grid(True)
    # plt.subplot(143)
    # plt.plot(uvw[:,2],-bins,'k',marker='o')
    # plt.plot(tuvw[:,2],-bins,'k',marker='+')
    # plt.grid(True)
    # plt.subplot(144)
    # plt.plot(uvw[:,3],-bins,'g',marker='o')
    # plt.plot(tuvw[:,3],-bins,'g',marker='+')
    # plt.grid(True)
    # plt.show()
    
    # print("Paused")
    # input()
    return tuvw

def qaqc_data(uvw, tEI, tC, tPG):
    #print('PROCESSING DVL DATA')
    corr_cut = 50
    ei_cut   = 70
    pg_cut   = 80
    
        
    # Change filled values to NaN
    uvw[uvw == -32768] = float("NAN")
    
    # Convert from mm/s to m/s
    uvw = uvw/1000

    uvw[:,0][tC[:,0] < corr_cut] = float("NAN")
    uvw[:,1][tC[:,1] < corr_cut] = float("NAN")
    uvw[:,2][tC[:,2] < corr_cut] = float("NAN")
    uvw[:,3][tC[:,3] < corr_cut] = float("NAN")

    uvw[:,0][tEI[:,0] < ei_cut] = float("NAN")
    uvw[:,1][tEI[:,1] < ei_cut] = float("NAN")
    uvw[:,2][tEI[:,2] < ei_cut] = float("NAN")
    uvw[:,3][tEI[:,3] < ei_cut] = float("NAN")
    
    uvw[:,0][tPG[:,0] < pg_cut] = float("NAN")
    uvw[:,1][tPG[:,1] < pg_cut] = float("NAN")
    uvw[:,2][tPG[:,2] < pg_cut] = float("NAN")
    uvw[:,3][tPG[:,3] < pg_cut] = float("NAN")


    # u1[c1 < corr_cut] = float("NAN")
    # u2[c2 < corr_cut] = float("NAN")
    # u3[c3 < corr_cut] = float("NAN")
    # u4[c4 < corr_cut] = float("NAN")  

    # u1[ei1 < ei_cut] = float("NAN")
    # u2[ei2 < ei_cut] = float("NAN")
    # u3[ei3 < ei_cut] = float("NAN")
    # u4[ei4 < ei_cut] = float("NAN")
    
    # u1[pg1 < pg_cut] = float("NAN")
    # u2[pg2 < pg_cut] = float("NAN")
    # u3[pg3 < pg_cut] = float("NAN")
    # u4[pg4 < pg_cut] = float("NAN")  
    
    return(uvw)

                        
# def process_data(U,V,H,dz,u_daverage,v_daverage):
#     global O_ls, G_ls, bin_new    
    
#     ## Feb-2021 jgradone@marine.rutgers.edu Initial
#     ## Jul-2021 jgradone@marine.rutgers.edu Updates for constraints
    
#     ## Purpose: Take velocity measurements from glider mounted ADCP and compute
#     # shear profiles
    
#     ## Outputs:
#     # O_ls is the ocean velocity profile
#     # G_ls is the glider velocity profile
#     # bin_new are the bin centers for the point in the profiles
#     # C is the constant used in the constraint equation (Not applicable for
#     # real-time processing)
    
#     ## Inputs:
#     # dz is desired vertical resolution, should not be smaller than bin length
#     # H is the max depth of the water column
#     # U is measured east-west velocities from ADCP
#     # V is measured north-south velocities from ADCP
#     # Z is the measurement depths of U and V
#     # uv_daverage is depth averaged velocity (Set to 0 for real-time)
    
#     ##########################################################################        
#     # Take difference between bin lengths for bin size [m]
#     bin_size = np.diff(bins)[0]
#     bin_num = len(bins)

#     # This creates a grid of the ACTUAL depths of the ADCP bins by adding the
#     # depths of the ADCP bins to the actual depth of the instrument
#     [bdepth,bbins]=np.meshgrid(depth,bins)
#     bin_depth = bdepth+bbins  
#     Z = bin_depth

#     # Set knowns from Equations 19 from Visbeck (2002) page 800
#     # Maximum number of observations (nd) is given by the number of velocity
#     # estimates per ping (nbin) times the number of profiles per cast (nt)
#     nbin = U.shape[0]  # number of programmed ADCP bins per individual profile
#     nt   = U.shape[1]  # number of individual velocity profiles
#     nd   = nbin*nt      # G dimension (1) 

#     # Define the edges of the bins
#     bin_edges = np.arange(0,math.floor(np.max(bin_depth)),dz).tolist()

#     # Check that each bin has data in it
#     bin_count = np.empty(len(bin_edges)-1) # Preallocate memory
#     bin_count[:] = np.NaN

#     for k in np.arange(len(bin_edges))[:-1]:
#         # Create index of depth values that fall inside the bin edges
#         ii = np.where((bin_depth > bin_edges[k]) & (bin_depth < bin_edges[k+1]))
#         bin_count[k] = len(bin_depth[ii])
#         ii = []

#     # Create list of bin centers    
#     bin_new = [x+dz/2 for x in bin_edges[:-1]]

#     # Chop off the top of profile if no data
#     ind = np.argmax(bin_count > 0) # Stops at first index greater than 0
#     bin_new = bin_new[ind:]        # Removes all bins above first with data
#     z1 = bin_new[0]                # Depth of center of first bin with data

#     # Create and populate G
#     nz = len(bin_new)  # number of ocean velocities desired in output profile
#     nm = nz + nt       # G dimension (2), number of unknowns
#     # Let's build the corresponding coefficient matrix G 
#     G = np.zeros((nd,nm))

#     # Indexing of the G matrix was taken from Todd et al. 2012
#     for ii in np.arange(nt):           # Number of ADCP profiles per cast
#         for jj in np.arange(nbin):     # Number of measured bins per profile

#             # Uctd part of matrix
#             G[(nbin*(ii-1))+jj,ii] = 1

#             # This will fill in the Uocean part of the matrix. It loops through
#             # all Z members and places them in the proper location in the G matrix

#             # Find the difference between all bin centers and the current Z value        
#             dx = abs(bin_new-Z[ii,jj])

#             # Find the minimum of these differences
#             minx = np.nanmin(dx)

#             # Finds bin_new index of the first match of Z and bin_new    
#             idx = np.argmin(dx-minx)

#             G[(nbin*(ii-1))+jj,nt+idx] = 1
#             del dx, minx, idx


#     # Reshape U and V into the format of the d column vector
#     d_u = np.flip(U.transpose(),axis=0)
#     d_u = d_u.flatten()
#     d_v = np.flip(V.transpose(),axis=0)
#     d_v = d_v.flatten()


#     ##########################################################################
#     ## This chunk of code containts the constraints for depth averaged currents
#     ## which we likely won't be using for the real-time processing

#     # Need to calculate C (Todd et al. 2017) based on our inputs 
#     # This creates a row that has the same # of columns as G. The elements
#     # of the row follow the trapezoid rule which is used because of the
#     # extension of the first bin with data to the surface. The last entry of
#     # the row corresponds to the max depth reached by the glider, any bins
#     # below that should have already been removed.

#     constraint = np.concatenate(([np.zeros(nt)], [z1/2], [z1/2+dz/2], [[dz]*(nz-3)], [dz/2]), axis=None)

#     # To find C, we use the equation of the norm and set norm=1 because we
#     # desire unity. The equation requires we take the sum of the squares of the
#     # entries in constraint.

#     sqr_constraint = constraint*constraint
#     sum_sqr_constraint = np.sum(sqr_constraint)

#     # Then we can solve for the value of C needed to maintain unity 

#     C = H*(1/np.sqrt(sum_sqr_constraint))

#     # This is where you would add the constraint for the depth averaged
#     # velocity from Todd et al., (2011/2017)

#     # OG
#     du = np.concatenate(([d_u],[C*u_daverage]), axis=None)
#     dv = np.concatenate(([d_v],[C*v_daverage]), axis=None)

#     # Build Gstar
#     # Keep this out because not using depth averaged currents
#     Gstar = np.vstack((G, (C/H)*constraint))

#     ##########################################################################

#     # Build the d matrix
#     d = list(map(complex,du, dv))

#     ##### Inversion!
#     ## If want to do with a sparse matrix sol'n, look at scipy
#     #Gs = scipy.sparse(Gstar)
#     Gs = Gstar

#     ms = np.linalg.solve(np.dot(Gs.conj().transpose(),Gs),Gs.conj().transpose())

#     ## This is a little clunky but I think the dot product fails because of
#     ## NaN's in the d vector. So, this code will replace NaN's with 0's just
#     ## for that calculation    
#     sol = np.dot(ms,np.where(np.isnan(d),0,d))

#     O_ls = sol[nt:]   # Ocean velocity
#     G_ls = sol[0:nt]  # Glider velocity





##################################################################################################

def inversion(U,V,dz,u_daverage,v_daverage,bins,depth, wDAC, wSmoothness):
    global O_ls, G_ls, bin_new

    ## Feb-2021 jgradone@marine.rutgers.edu Initial
    ## Jul-2021 jgradone@marine.rutgers.edu Updates for constraints
    ## Jun-2022 jgradone@marine.rutgers.edu Corrected dimensions and indexing of G matrix
    ## Jun-2022 jgradone@marine.rutgers.edu Added curvature minimizing constraint and constraint weights

    ## Purpose: Take velocity measurements from glider mounted ADCP and compute
    # shear profiles

    ## Outputs:
    # O_ls is the ocean velocity profile
    # G_ls is the glider velocity profile
    # bin_new are the bin centers for the point in the profiles
    # obs_per_bin is the number of good velocity observations per final profile bin

    ## Inputs:
    # dz is desired vertical resolution, should not be smaller than bin length 
    # U is measured east-west velocities from ADCP
    # V is measured north-south velocities from ADCP
    # bins is the bin depths for the U and V measurements
    # uv_daverage is depth averaged velocity (Set to 0 for real-time)
    # depth is the depth of the glider measured by the ADCP
    # wDAC is the weight of the DAC constraint (5 per Todd et al. 2017)
    # wSmoothness is the weight of the curvature minimizing contraint (1 per Todd et al. 2017)


    #########################################################################  
    ## These steps filter for NAN rows and columns so they are technically QAQC
    ## but I think the best place to put them is inthe inversion function because
    ## if there are NaNs still present in the data here, it will throw everything off
    ## These steps are HUGE for efficiency because it reduces the size of the G
    ## matrix as much as possible.

    ## This determines the rows (bins) where all the columns are NaN
    nanind = np.where( (np.sum(np.isnan(U),axis=1)/U.shape[1]) == 1)[0]
    if len(nanind) > 0:
        U = np.delete(U,nanind,axis=0)
        V = np.delete(V,nanind,axis=0)
        bins = np.delete(bins,nanind)

    ## Do the same thing with individual ensembles. Note: need to remove the corresponding
    ## ensemble pressure reading to ensure correction dimensions and values.
    nanind = np.where((np.sum(np.isnan(U),axis=0)/U.shape[0]) == 1)[0]
    if len(nanind) > 0:
        U = np.delete(U,nanind,axis=1)
        V = np.delete(V,nanind,axis=1)
        depth = np.delete(depth,nanind)
    ##########################################################################        


    ##########################################################################        
    # Take difference between bin lengths for bin size [m]
    bin_size = np.diff(bins)[0]
    bin_num = len(bins)
    # This creates a grid of the ACTUAL depths of the ADCP bins by adding the
    # depths of the ADCP bins to the actual depth of the instrument
    [bdepth,bbins]=np.meshgrid(depth,bins)
    bin_depth = bdepth+bbins  
    Z = bin_depth
    # Calculate the maximum depth of glider which is different than maximum ADCP bin depth
    ZmM = np.nanmax(depth)
    ##########################################################################        


    ##########################################################################        
    # Set knowns from Equations 19 from Visbeck (2002) page 800
    # Maximum number of observations (nd) is given by the number of velocity
    # estimates per ping (nbin) times the number of profiles per cast (nt)
    nbin = U.shape[0]  # number of programmed ADCP bins per individual profile
    nt   = U.shape[1]  # number of individual velocity profiles
    nd   = nbin*nt      # G dimension (1) 

    # Define the edges of the bins
    bin_edges = np.arange(0,math.floor(np.max(bin_depth)),dz).tolist()

    # Check that each bin has data in it
    bin_count = np.empty(len(bin_edges)-1) # Preallocate memory
    bin_count[:] = np.NaN

    for k in np.arange(len(bin_edges))[:-1]:
        # Create index of depth values that fall inside the bin edges
        ii = np.where((bin_depth > bin_edges[k]) & (bin_depth < bin_edges[k+1]))
        bin_count[k] = len(bin_depth[ii])
        ii = []

    # Create list of bin centers    
    bin_new = [x+dz/2 for x in bin_edges[:-1]]

    # Calculate which FINAL solution bin is deeper than the maximum depth of the glider
    # This is done so that the depth averaged velocity constraint is only applied to bins shallower than this depth
    depth_ind = len(np.where(bin_new>ZmM)[0])
    # Chop off the top of profile if no data
    ind = np.argmax(bin_count > 0) # Stops at first index greater than 0
    bin_new = bin_new[ind:]        # Removes all bins above first with data
    z1 = bin_new[0]                # Depth of center of first bin with data
    ##########################################################################        


    ##########################################################################        
    # Create and populate G
    nz = len(bin_new)  # number of ocean velocities desired in output profile
    nm = nt + nz       # G dimension (2), number of unknowns
    # Let's build the corresponding coefficient matrix G 
    G = scipy.sparse.lil_matrix((nd, nm), dtype=float)

    # Indexing of the G matrix was taken from Todd et al. 2012
    for ii in np.arange(0,nt):           # Number of ADCP ensembles per segment
        for jj in np.arange(0,nbin):     # Number of measured bins per ensemble 

            # Uctd part of matrix
            G[(nbin*(ii))+jj,ii] = -1
            # This will fill in the Uocean part of the matrix. It loops through
            # all Z members and places them in the proper location in the G matrix
            # Find the difference between all bin centers and the current Z value        
            dx = abs(bin_new-Z[jj,ii])
            # Find the minimum of these differences
            minx = np.nanmin(dx)
            # Finds bin_new index of the first match of Z and bin_new    
            idx = np.argmin(dx-minx)

            # Uocean part of matrix
            G[(nbin*(ii))+jj,(nt)+idx] = 1

            del dx, minx, idx

    ##########################################################################        
    # Reshape U and V into the format of the d column vector (order='F')
    # Based on how G is made, d needs to be ensembles stacked on one another vertically
    d_u = U.flatten(order='F')
    d_v = V.flatten(order='F')

    ##########################################################################
    ## This chunk of code containts the constraints for depth averaged currents
    # Make sure the constraint is only applied to the final ocean velocity bins that the glider dives through
    # Don't apply it to the first bin and don't apply it to the bins below the gliders dive depth
    constraint = np.concatenate(([np.zeros(nt)], [0], [np.tile(dz,nz-(1+depth_ind))], [np.zeros(depth_ind)]), axis=None)

    # Ensure the L^2 norm of the constraint equation is unity
    constraint_norm = np.linalg.norm(constraint/ZmM)
    C = 1/constraint_norm
    constraint_normalized = (C/ZmM)*constraint ## This is now equal to 1 (unity)
    # Build Gstar and add weight from todd 2017
    ## Some smarts would be to calculate signal to noise ratio first
    Gstar = scipy.sparse.vstack((G,wDAC*constraint_normalized), dtype=float)


    # Add the constraint for the depth averaged velocity from Todd et al. (2017)
    du = np.concatenate(([d_u],[wDAC*C*u_daverage]), axis=None)
    dv = np.concatenate(([d_v],[wDAC*C*v_daverage]), axis=None)
    d = np.array(list(map(complex,du, dv)))


    ##########################################################################        
    #### THIS removes all NaN elements of d AND Gstar so the inversion doesn't blow up with NaNs
    ind2 = np.where(np.isnan(d)==True)[0]
    d = np.delete(d,ind2)

    def delete_rows_csr(mat, indices):
        """
        Remove the rows denoted by ``indices`` form the CSR sparse matrix ``mat``.
        """
        if not isinstance(mat, scipy.sparse.csr_matrix):
            raise ValueError("works only for CSR format -- use .tocsr() first")
        indices = list(indices)
        mask = np.ones(mat.shape[0], dtype=bool)
        mask[indices] = False
        return mat[mask]

    Gstar = delete_rows_csr(Gstar.tocsr().copy(),ind2)

    #########################################################################        
    # Test adding depth for tracking bin location
    # d is ensembles stacked on one another vertically so same for Z (order='F')
    Z_filt = Z.flatten(order='F')
    Z_filt = np.delete(Z_filt,ind2)
    Z_filt = np.concatenate(([Z_filt],[0]), axis=None)

    ##########################################################################        
    ## Calculation the number of observations per bin
    obs_per_bin = np.empty(len(bin_new))
    obs_per_bin[:] = np.NaN

    for x in np.arange(0,nz):
        rows_where_nt_not_equal_zero = np.where(Gstar.tocsr()[0:Z_filt.shape[0],nt+x].toarray() > 0 )[0]
        obs_per_bin[x] = len(rows_where_nt_not_equal_zero)

    ## If there is no data in the last bin, drop that from the G matrix, bin_new, and obs_per_bin
    if obs_per_bin[-1] == 0:
        Gstar.tocsr()[:,:-1]
        bin_new = bin_new[:-1]
        obs_per_bin = obs_per_bin[:-1]
        ## Update nz and nt
        nz = len(bin_new)
        nt = Gstar.shape[1]-nz

    ##########################################################################        
    ## Smoothness constraint
    ## Only do this is the smoothness constraint is set
    if wSmoothness > 0:
        ## Add a vector of zerosm the length of nz, twice to the bottom of the data column vector
        d = np.concatenate(([d],[np.zeros(nz)],[np.zeros(nz)]), axis=None)
        ## Constraint on smoothing Uocean side of matrix
        smoothing_matrix_Uocean = scipy.sparse.diags([[-1],[2],[-1]], [0,1,2], shape=(nz,nz))
        smoothing_matrix1 = scipy.sparse.hstack((np.zeros((nz,nt)),smoothing_matrix_Uocean), dtype=float)
        ## Constraint on smoothing Uglider side of matrix
        smoothing_matrix_Uglider = scipy.sparse.diags([[-1],[2],[-1]], [0,1,2], shape=(nz,nt))
        smoothing_matrix2 = scipy.sparse.hstack((smoothing_matrix_Uglider,np.zeros((nz,nz))), dtype=float)
        Gstar = scipy.sparse.vstack((Gstar,wSmoothness*smoothing_matrix1,wSmoothness*smoothing_matrix2), dtype=float)


    ##########################################################################        
    ## Run the Least-Squares Inversion!
    x = lsqr(Gstar, d)[0]

    O_ls = x[nt:]
    G_ls = x[0:nt] 
    ########################################################################## 

    return(O_ls, G_ls, bin_new, obs_per_bin)





    
    
    
    
    
def write_data(infile):
    global O_ls, G_ls, bin_new,time
    basefile=os.path.basename(infile).upper()
    ncfile=odir+basefile.replace('PD0','nc')
    print('Writing profile DATA : '+ncfile)
    ncfile = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')
    ncfile.Conventions= "CF-1.6"
#    
#    
    ncfile.createDimension('time', 1)
    ncfile.createDimension('depth', len(bin_new))
    ncfile.createVariable('time','f8',('time'))
    ncfile.createVariable('timemax','f8',('time'))
    ncfile.createVariable('timemin','f8',('time'))
    ncfile.createVariable('depth','f8',('depth'))
    
    ncfile.createVariable('u','f8',('time','depth'))
    ncfile.createVariable('v','f8',('time','depth'))
    print(np.mean(time))
    
# #        

    ncfile.variables['time'][:]=np.mean(time)
    ncfile.variables['time'].units='days since 2020-01-01 00:00:00'
    ncfile.variables['time'].long_name='time'
    ncfile.variables['time'].standard_name='time'


    ncfile.variables['timemax'][:]=np.max(time)
    ncfile.variables['timemax'].units='days since 2020-01-01 00:00:00'
    ncfile.variables['timemax'].long_name='maximum time'
    ncfile.variables['timemax'].standard_name='timemax'


    ncfile.variables['timemin'][:]=np.min(time)
    ncfile.variables['timemin'].units='days since 2020-01-01 00:00:00'
    ncfile.variables['timemin'].long_name='minimum time'
    ncfile.variables['timemin'].standard_name='timemin'

    
# #    
    ncfile.variables['depth'][:]=bin_new
    ncfile.variables['depth'].units='m'
    ncfile.variables['depth'].long_name='Cell Depth'
    ncfile.variables['depth'].standard_name='depth'
#        plt.plot(np.real(O_ls),bin_new,label='u - velocity')
#    plt.plot(np.imag(O_ls),bin_new,label='v - velocity')

    ncfile.variables['u'][:]=np.real(O_ls)
    ncfile.variables['u'].units='m s-1'
    ncfile.variables['u'].long_name='eastward water velocity'
    ncfile.variables['u'].standard_name='eastward_sea_water_velocity'

    ncfile.variables['v'][:]=np.imag(O_ls)
    ncfile.variables['v'].units='m s-1'
    ncfile.variables['v'].long_name='northward water velocity'
    ncfile.variables['v'].standard_name='northward_sea_water_velocity'

    # ncfile.variables['w'][:]=w/1000.0
    # ncfile.variables['w'].units='m s-1'
    # ncfile.variables['w'].long_name='vertical water velocity'
    # ncfile.variables['w'].standard_name='upward_sea_water_velocity'

    # ncfile.variables['err'][:]=err/1000.0
    # ncfile.variables['err'].units='m s-1'
    # ncfile.variables['err'].long_name='error water velocity'
    
    # print u.shape
    # print ei1.shape
    # ncfile.variables['ei1'][:]=ei1
    # ncfile.variables['ei1'].units='counts'
    # ncfile.variables['ei1'].long_name='beam 1 echo intensity'
    
    # ncfile.variables['ei2'][:]=ei2
    # ncfile.variables['ei2'].units='counts'
    # ncfile.variables['ei2'].long_name='beam 1 echo intensity'
    
    # ncfile.variables['ei3'][:]=ei3
    # ncfile.variables['ei3'].units='counts'
    # ncfile.variables['ei3'].long_name='beam 1 echo intensity'
    
    # ncfile.variables['ei4'][:]=ei4
    # ncfile.variables['ei4'].units='counts'
    # ncfile.variables['ei4'].long_name='beam 1 echo intensity'
            
            
    
#    
    ncfile.close()
def sph2cart(az,elev,r):
    z = r * np.sin(elev)
    rcoselev = r * np.cos(elev)
    x = rcoselev * np.cos(az)
    y = rcoselev * np.sin(az)
    return x,y,z;


def plot_data():
    print('Plotting DVL DATA')       
    # print(bins)
    # plt.figure(1)
    # plt.clf()
    # plt.plot(time,-depth,'r')
    # plt.ylabel('Depth [m]')
    # plt.grid(True)
    
    plt.figure(2)
    plt.clf()
    plt.subplot(311)
    plt.plot(time,heading,'r')
    plt.ylabel('Heading')
    plt.grid(True)
    plt.subplot(312)
    plt.plot(time,pitch,'r')
    plt.ylabel('Pitch')
    plt.grid(True)
    plt.subplot(313)
    plt.plot(time,roll,'r')
    plt.ylabel('Roll')
    plt.grid(True)
    
    
    cmap = plt.get_cmap('jet')
    [x,y]=np.meshgrid(time,bins)
    [bdepth,bbins]=np.meshgrid(depth,bins)

    by=bdepth+bbins
    

    
    
    
    fig3=plt.figure(5)
    plt.clf()
    
    ax1=plt.subplot(411)
    pc2=plt.pcolormesh(time,-bins,u1.transpose(),cmap=cmap,vmin=-1,vmax=1)
    #plt.plot(time,-depth,'k')
    fig3.colorbar(pc2,ax=ax1)

    ax1=plt.subplot(412)
    pc2=plt.pcolormesh(time,-bins,u2.transpose(),cmap=cmap,vmin=-1,vmax=1)
    #plt.plot(time,-depth,'k')
    fig3.colorbar(pc2,ax=ax1)
    
    ax1=plt.subplot(413)
    pc2=plt.pcolormesh(time,-bins,u3.transpose(),cmap=cmap,vmin=-1,vmax=1)
    #plt.plot(time,-depth,'k')
    fig3.colorbar(pc2,ax=ax1)
    
    ax1=plt.subplot(414)
    pc2=plt.pcolormesh(time,-bins,u4.transpose(),cmap=cmap,vmin=-1,vmax=1)
    #plt.plot(time,-depth,'k')
    fig3.colorbar(pc2,ax=ax1)


    fig4=plt.figure(6)
    plt.clf()
    ax1=plt.subplot(411)
    pc2=plt.pcolormesh(time,-bins,c1.transpose(),cmap=cmap,vmin=0,vmax=80)
    #plt.plot(time,-depth,'k')
    fig4.colorbar(pc2,ax=ax1)

    ax1=plt.subplot(412)
    pc2=plt.pcolormesh(time,-bins,c2.transpose(),cmap=cmap,vmin=0,vmax=80)
    #plt.plot(time,-depth,'k')
    fig4.colorbar(pc2,ax=ax1)
   
    ax1=plt.subplot(413)
    pc2=plt.pcolormesh(time,-bins,c3.transpose(),cmap=cmap,vmin=0,vmax=80)
    #plt.plot(time,-depth,'k')
    fig4.colorbar(pc2,ax=ax1)
   
    ax1=plt.subplot(414)
    pc2=plt.pcolormesh(time,-bins,c4.transpose(),cmap=cmap,vmin=0,vmax=80)
    #plt.plot(time,-depth,'k')
    fig4.colorbar(pc2,ax=ax1)
    ## A few test plots
    
    

    fig4=plt.figure(7)
    plt.clf()
    ax1=plt.subplot(411)
    pc2=plt.pcolormesh(time,-bins,ei1.transpose(),cmap=cmap,vmin=0,vmax=80)
    #plt.plot(time,-depth,'k')
    fig4.colorbar(pc2,ax=ax1)

    ax1=plt.subplot(412)
    pc2=plt.pcolormesh(time,-bins,ei2.transpose(),cmap=cmap,vmin=0,vmax=80)
    #plt.plot(time,-depth,'k')
    fig4.colorbar(pc2,ax=ax1)
   
    ax1=plt.subplot(413)
    pc2=plt.pcolormesh(time,-bins,ei3.transpose(),cmap=cmap,vmin=0,vmax=80)
    #plt.plot(time,-depth,'k')
    fig4.colorbar(pc2,ax=ax1)
   
    ax1=plt.subplot(414)
    pc2=plt.pcolormesh(time,-bins,ei4.transpose(),cmap=cmap,vmin=0,vmax=80)
    #plt.plot(time,-depth,'k')
    fig4.colorbar(pc2,ax=ax1)
    ## A few test plots
    
    
    plt.figure(8)
    plt.plot(np.real(O_ls),bin_new,label='u - velocity')
    plt.plot(np.imag(O_ls),bin_new,label='v - velocity')
    plt.ylim(1000,0)
    plt.legend()
    # plt.figure(9)
    # plt.plot(np.real(G_ls),bin_new,label='u - velocity')
    # plt.plot(np.imag(G_ls),bin_new,label='v - velocity')
    # plt.ylim(30,0)
    # plt.legend()
  
    # # Just give me one profile
    # fig3 = plt.figure(7)
    # pc3 = plt.pcolormesh(x,-by,u1.transpose(),cmap=cmap,vmin=-1,vmax=1)
    # plt.plot(time,-depth,'k')
    # plt.title('20 < abs(Pitch) < 30')    
    # fig3.colorbar(pc3)
    
    # dp = np.diff(pitch)
    # plt.figure(8)
    # plt.plot(dp[1:200])
    # plt.axhline(y=-5)
    
#def bin(s):
#    return str(s) if s<=1 else bin(s>>1) + str(s&1)
    
if __name__ == "__main__":
    

    
    print('RUNNING')
    main(sys.argv)
    print('FINISHED')
    

