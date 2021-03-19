#!/usr/bin/env python
#
import os,sys
import glob
import numpy as np
import matplotlib,netCDF4
import matplotlib.pyplot as plt
import datetime
import struct
#
rtime=datetime.datetime(2013,1,1,0,0,0)

odir='./'
idir='C:\\work\\glideradcp\\data\\ru33_2020_11_20_dvl\\pd0\\'

def main(argv):
    
    files=glob.glob(idir+'*.PD0')
    files.sort(key=os.path.getmtime)
    file=files[-1]
    print(file)
    read_write_PD0(file)
    
def read_write_PD0(infile):
    print('PROCESSING_ADCP')
    #print(infile)


    f=open(infile,'rb')
    dat = f.read()
    f.close()

    for [ind,tmp1] in enumerate(dat):
          if hex(dat[ind])=='0x7f' and hex(dat[ind+1]) =='0x7f':
                 
                  break

    nbytes=struct.unpack("h",dat[ind+2:ind+4])[0]+2
    
    
# ##    
    print('Finding Ensembles')      
  #  print(nbytes)
    Iens=[]
    nind=0
    n=0
    for [ind,tmp1] in enumerate(dat):
          if hex(dat[ind])=='0x7f' and hex(dat[ind+1]) =='0x7f':
           
              n=n+1
         #    a=buffer(dat,ind+2,2)

              nbytes2=struct.unpack("h",dat[ind+2:ind+4])[0]+2  
# #             
             
              startens=ind
              tdat=dat[startens:startens+nbytes]
              if len(tdat)<nbytes:
                   print('breaking')
                   break
#               #a=buffer(tdat,nbytes-2,2)
              
              tmp=tdat[nbytes-2:nbytes]
              chksum=struct.unpack("<H",tmp)[0]
              if (sum(tdat[:nbytes-2]) & 65535) ==  chksum:
                      
                   if nbytes == nbytes2:
  
                       nind=ind
                       Iens.append(ind)
 #             else:
 #                print('Bad Checksum')
#                 
 #   print('+++++')   

 #   print(len(Iens))
                 
# #             
    nens=len(Iens) 
#   print(nens)
    print('creating arrays')
    time=np.empty((nens),np.float)    
    pressure=np.empty((nens),np.float) 
    pitch=np.empty((nens),np.float) 
    roll=np.empty((nens),np.float) 
    heading=np.empty((nens),np.float) 
    temp=np.empty((nens),np.float) 
#        
    u=np.empty((nens,100),np.float) 
    v=np.empty((nens,100),np.float) 
    w=np.empty((nens,100),np.float)
    err=np.empty((nens,100),np.float)  
    ei1=np.empty((nens,100),np.float) 
    ei2=np.empty((nens,100),np.float) 
    ei3=np.empty((nens,100),np.float)
    ei4=np.empty((nens,100),np.float)  
    c1=np.empty((nens,100),np.float) 
    c2=np.empty((nens,100),np.float) 
    c3=np.empty((nens,100),np.float)
    c4=np.empty((nens,100),np.float)   
    
    #ADDDED 02/11/2020
    pg1=np.empty((nens,100),np.float) 
    pg2=np.empty((nens,100),np.float) 
    pg3=np.empty((nens,100),np.float)
    pg4=np.empty((nens,100),np.float)              
 #   print('loop')
    ind=0
    eoffset=0
    for ind2 in Iens:
#        print '---'
        startens=(ind2)
        tdat=dat[startens:startens+nbytes]
        # a=buffer(tdat,2,2)
        tnbytes=struct.unpack("H",tdat[2:4])[0]+2
        # a=buffer(tdat,nbytes-2,2)
        chksum=struct.unpack("<H",tdat[nbytes-2:nbytes])[0]

        if (sum(tdat[:nbytes-2]) & 65535) ==  chksum:
 #             print("stuff")
              ndtype=struct.unpack("b",tdat[5:6])[0]
 
              offsets=list()
              for ind3 in range(ndtype):
                   Is=6+ind3*2
                   offsets.append(struct.unpack_from("h",tdat[Is:Is+2])[0])
                
                     
#             #FIXEDLEADER
# #            print ind2
              Is=offsets[0]+8
              nbeam=tdat[Is]
                
              Is=offsets[0]+9
              ncells=tdat[Is]  

              Is=offsets[0]+12
              cellsize=struct.unpack("H",tdat[Is:Is+2])[0]        
              cellsize=cellsize/100.0       

              Is=offsets[0]+32
              bin1=struct.unpack("H",tdat[Is:Is+2])[0]    
              bin1=bin1/100.0
     #         print(cellsize)
      #        print(bin1)
              
              Is=offsets[1]+2
              ens=struct.unpack("H",tdat[Is:Is+2])[0]        
        #      print(ens)
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
          #   print [year+2000,month,day,hour,minute,sec,hsec*10]
              ttime = datetime.datetime(year+2000,month,day,hour,minute,sec,hsec*10)-rtime
#            trtime = datetime.datetime(year+2000,month,day,hour,minute,sec,hsec)
     #         print([year+2000,month,day,hour,minute,sec])
              Is=offsets[1]+18
              theading=struct.unpack("H",tdat[Is:Is+2])[0]        
              theading=theading/100.0    
              Is=offsets[1]+20
              tpitch=struct.unpack("h",tdat[Is:Is+2])[0]        
              tpitch=tpitch/100.0    
              Is=offsets[1]+22
              troll=struct.unpack("h",tdat[Is:Is+2])[0]        
              troll=troll/100.0    
            
              Is=offsets[1]+26
              ttemp=struct.unpack("h",tdat[Is:Is+2])[0]        
              ttemp=ttemp/100.0       
            
              Is=offsets[1]+48
              tpress=struct.unpack("i",tdat[Is:Is+4])[0]        
              tpress=tpress/1000.0       
            
            
   #           print([theading,tpitch,troll,ttemp,tpress])
    #          print(tpress)
#             #velocity data 
#             a=buffer(tdat,offsets[2]+2,ncells*4*2)
# #            print ncells*4*2
              Is=offsets[2]+2
              fmt = "<%dh" % (ncells*4)
              uvw=struct.unpack(fmt,tdat[Is:Is+ncells*4*2])
              uvw=np.array(uvw)

#             #EI data 

              Is=offsets[3]+2
              fmt = "<%dB" % (ncells*4)
              tEI=struct.unpack(fmt,tdat[Is:Is+ncells*4])
              tEI=np.array(tEI)
        

#             #C data 
              Is=offsets[4]+2
              fmt = "<%dB" % (ncells*4)
              tC=struct.unpack(fmt,tdat[Is:Is+ncells*4])
              tC=np.array(tC)
              
              #ADDDED 02/11/2020
#             #PG data 
              Is=offsets[5]+2
              fmt = "<%dB" % (ncells*4)
              tPG=struct.unpack(fmt,tdat[Is:Is+ncells*4])
              tPG=np.array(tPG)
              
            
              uvw.shape=(ncells,4)
              tEI.shape=(ncells,4)
              tC.shape=(ncells,4)
              tPG.shape=(ncells,4)#ADDDED 02/11/2020
             
             
             
    #          print(uvw)

# #            print ttime.days
              time[ind]=ttime.days+ttime.seconds/86400.0
      
              pressure[ind]=tpress
              pitch[ind]=tpitch
              roll[ind]=troll
              temp[ind]=ttemp
              heading[ind]=theading
              u[ind,0:ncells]=uvw[:,0]
              v[ind,0:ncells]=uvw[:,1]
              w[ind,0:ncells]=uvw[:,2]
              err[ind,0:ncells]=uvw[:,3]
        
              ei1[ind,0:ncells]=tEI[:,0]
              ei2[ind,0:ncells]=tEI[:,1]
              ei3[ind,0:ncells]=tEI[:,2]
              ei4[ind,0:ncells]=tEI[:,3]
        
              c1[ind,0:ncells]=tC[:,0]
              c2[ind,0:ncells]=tC[:,1]
              c3[ind,0:ncells]=tC[:,2]
              c4[ind,0:ncells]=tC[:,3]
              
              #ADDDED 02/11/2020
              pg1[ind,0:ncells]=tPG[:,0]
              pg2[ind,0:ncells]=tPG[:,1]
              pg3[ind,0:ncells]=tPG[:,2]
              pg4[ind,0:ncells]=tPG[:,3]
              
              
              
              
              ind=ind+1
        else:
          #   print 'BAD CHECKSUM'
            # eoffset=eoffset+1
            
            continue
     
    
# #            
#    nt=len(time)
#    print(nt)


    
    print('FINISHED')
        
    
    
if __name__ == "__main__":
    

    
    print('RUNNING')
    main(sys.argv)
