#!/usr/bin/env python
#upgrading to scipy version 0.100dev , changing interface
import numpy as np
import numpy.fft as ft
import scipy.fftpack as sft
from scipy.io import netcdf_file
import matplotlib.pyplot as plt
import os
from matplotlib import ticker, cm
import fortranformat as ff
import f90nml

#other deps below
from periodictable import elements

def print_vector(nrep,fstr,a):
    """
    Converts an array of numbers into a string formated by fstr with
    nrep values per line.
    """
    n=a.size
    pa=""
    ta=a.reshape((a.size),order='F')
    for j in range(0,n,nrep):
        pa=pa+ "".join(map(lambda f: fstr % f, ta[j:min(j+nrep,n)]))+"\n"
    return pa


def ListToFormattedString(alist,fstr):
    # Create a format spec for each item in the input `alist`.
    # E.g., each item will be right-adjusted, field width=3.
    format_list = [fstr for item in alist]

    # Now join the format specs into a single string:
    # E.g., '{:>3}, {:>3}, {:>3}' if the input list has 3 items.
    s = ', '.join(format_list)

    # Now unpack the input list `alist` into the format string. Done!
    return s.format(*alist)


def formattedwrite(file,a):
  sza=len(a)
  for idx in range(0,int(sza/4)*4,4):
    file.write(f"{a[idx]:18.9E}{a[idx+1]:18.9E}{a[idx+2]:18.9E}{a[idx+3]:18.9E}\n")
  rem = np.mod(sza,4)
  if rem>0: #print remaining elements
    file.write(''.join([ "%18.9E" % x for x in a[-rem:] ])+"\n") 


# Plasma species in template namelist
def get_spec_toric(toricnml):
    "Collect info on species for ICRF sim in TORIC in nice readable format"

    spec_toric=list(zip(map(round,toricnml['equidata']['atm']),
                        map(int,toricnml['equidata']['azi'])))

    #if idprof=1, read profiles and conc from text file instead.
    for i,s in enumerate(spec_toric):
      name=str(elements[s[1]][s[0]])
      if False:
          print('name',i,name,len(100*toricnml['equidata']['aconc']),
                100*toricnml['equidata']['aconc'][i] )
      spec_toric[i]={'name':name,'A':spec_toric[i][0],'Z':spec_toric[i][1],
                     'Conc%':100*toricnml['equidata']['aconc'][i]}
    spec_toric.insert(0,{'name':'e', 'A':0, 'Z':-1 , 'Conc%': 100})
    #print('spec',len(spec_toric),toricnml['equidata']['atm'],spec_toric)
    return spec_toric


def write_profnt(fname,equidt):
    """
    Inputs:
        fname: file to save
        equidt: python dictionary of profiles to write out
            rhopro: sqrt(Psipol/Psipol[a])
            ne: [cm-3] on rhopro mesh
            te: [keV]  on rhopro mesh
            iatm: array of atomic masses/ (C12/12)
            zi: array of atomic numbers
            ni: ion densities on rhopro mesh or concentration
            ti: ion temperatures on rhopro mesh

         version: Generally format2 is used. File is self describing in
                  number of elements and profiles.
    """

    import fortranformat as ff

    keys=['rhopro','ne','te']
    f0001=ff.FortranRecordWriter('a10,5i4')
    f0002=ff.FortranRecordWriter('2i4')
    f0003=ff.FortranRecordWriter('a10,1i4')
    f0004=ff.FortranRecordWriter('5e16.9')
    f0005=ff.FortranRecordWriter('a32')
    f0006=ff.FortranRecordWriter('1e16.9')
    
    with open(fname,'w') as of:
        iatm=equidt['iatm']
        iazi=equidt['iazi']
        nspec=equidt['nspec']
        mainsp=equidt['mainsp']
        aconc=equidt['aconc']
        nprodt=equidt['nprodt']
        kdiff_itemp=equidt['kdiff_itemp']
        kdiff_idens=equidt['kdiff_idens']
        variantid=equidt['variant']
        ion_temp=equidt['ion_temp']
        of.write( f0001.write([variantid,nprodt,nspec,mainsp,kdiff_idens,kdiff_itemp]) )
        if equidt['variant']=='Rfxqlo_Pro':
            rfxqlo_pro_variant = True
        else:
            rfxqlo_pro_variant = False

            
        for isp in range(nspec):
             of.write( f0002.write( [iatm[isp],iazi[isp]] ) )

        for profile in keys: #write electron profiles
            of.write( f0004.write( equidt[profile] ) )

             
            #mainsp=1
            #namelist['equidata']['mainsp']=mainsp
            kdiff_idens=equidt['kdiff_idens'] #0 #specify concentrations
            kdiff_itemp=equidt['kdiff_itemp'] #0 #one ion temp
            of.write('{:<10s}{:4d}{:4d}{:4d}{:4d}{:4d}\n'.
            format('profnt_py', nprodt,nspec, mainsp,kdiff_idens,kdiff_itemp))
            for isp in range(nspec):
                of.write('{:4d}{:4d}\n'.
                  format(int(equidt['iatm'][isp]),int(equidt['iazi'][isp])) )
            profiles=['rhopro','tbne','tbte']
            for profile in profiles:
                of.write('{:<10s}{:4d}\n'.format(profile, nprodt ))
                of.write(print_vector(5,'%16.9e',np.array(equidt[profile])))

            #write ion densities and temperatures
            for isp in range(nspec):
                if kdiff_idens==0:
                    of.write('{:<10s}\n'.format('ni_conc'+str(isp)))
                    of.write('%16.9e \n' % equidt['aconc'][isp])
                else:
                    of.write('{:<10s}\n'.format('tbni'+str(isp)))
                    of.write(print_vector(5,'%16.9e',equidt['ni'][:,isp]))

                if kdiff_itemp==0 and isp==0:
                    of.write('{:<10s}\n'.format('ion_temp') )
                    of.write(print_vector(5,'%16.9e',equidt['ti_provv']))

                if kdiff_itemp==1:
                    of.write('{:<10s}\n'.format('ion_temp'+str(isp)) )
                    of.write(print_vector(5,'%16.9e',equidt['ti_provv'][:,isp]))


def readArray(of,fmt,shp,nperline=5):
    vals=[]
    if len(shp)==1: N=shp[0]
    if len(shp)==2: N=shp[0]*shp[1]
    nlines = int(N/nperpline)
    if (N%nperpline)!=0: nlines+=1
    for i in range( nlines ):
        vals.extend(fmt.read(next(of)))
    return np.reshape(np.array(vals[0:N]),shp)
                    
                    
def read_equidt(filename,idebug=False):
    import fortranformat as ff

    keys=['rhopro','ne','te']
    profiles=dict.fromkeys(keys)
    f0001=ff.FortranRecordReader('a10,5i4')
    f0002=ff.FortranRecordReader('2i4')
    f0003=ff.FortranRecordReader('a10,1i4')
    f0004=ff.FortranRecordReader('5e16.9')
    f0005=ff.FortranRecordReader('a32')
    f0006=ff.FortranRecordReader('1e16.9')


    with open(filename,'r') as of:
        [variantid,nprodt,nspec,mainsp,kdiff_idens,kdiff_itemp]=f0001.read(next(of))
        if idebug: print('variantid',variantid,nprodt,nspec,mainsp,kdiff_idens,kdiff_itemp)
        if variantid=='Rfxqlo_Pro':
            #  Dmc -- define flag to indicate TRANSP format variant
            if idebug: print(' Detected: TRANSP "Rfxqlo_Pro" file format variant!')
            kdiff_itemp = 1
            kdiff_idens = 1
            rfxqlo_pro_variant = True
        else:
            rfxqlo_pro_variant = False

        if mainsp<=0 or mainsp>nspec:
            mainsp = 1

        if kdiff_itemp==0:
            nsptmp=1
        else:
            nsptmp=nspec

        iatm=np.zeros(nspec, dtype=int)
        iazi=np.zeros(nspec, dtype=int)
        aconc=np.zeros(nspec, dtype=float)    
        ion_temp=np.zeros([nprodt,nspec], dtype=float)
        ion_dens=np.zeros([nprodt,nspec], dtype=float)

        for isp in range(nspec):
            [iatm[isp],iazi[isp]] = f0002.read(next(of))
            #print('spec',int(iatm[isp]),int(iazi[isp]) )

        for profile in keys: #reading electron profiles
            [proname]=f0005.read(next(of))
            profiles[profile]=readArray(of,f0004,[nprodt])
            if idebug: print('reading name,size',proname,profiles[proname.strip()][0:4])

        if rfxqlo_pro_variant:  #Transp variant
            for isp in range(nspec):
                if isp==mainsp: continue
                [proname]=f0005.read(next(of))
                #print('reading ', proname)
                ion_dens[:,isp]=readArray(of,f0004,[nprodt])
            
            for isp in range(nspec):
                [proname]=f0005.read(next(of))
                if idebug: print('reading ', proname)
                ion_temp[:,isp]=readArray(of,f0004,[nprodt])

        else:
            for isp in range(nspec): # main TORIC convention
        
                if kdiff_idens != 0: # Different profiles for different species
                    [proname]=f0005.read(next(of))
                    ion_dens[:,isp]=readArray(of,f0004,[nprodt])

                if kdiff_idens < 0:
                    ion_dens[:,isp] *= profiles['ne']
                
                if kdiff_idens == 0: #scalar concentrations used
                    [proname]=f0005.read(next(of))
                    [aconc[isp]]=f0006.read(next(of))
                    if idebug: print('reading name',proname,isp,aconc[isp])

                if isp < nsptmp:
                    [proname]=f0005.read(next(of))
                    if idebug: print('reading ', proname)            
                    ion_temp[:,isp]=readArray(of,f0004,[nprodt])
    

    profiles['iatm'] =iatm    #atomic mass number
    profiles['iazi'] =iazi    #atomic charge number
    profiles['nspec'] =nspec
    profiles['mainsp'] =mainsp
    profiles['aconc'] =aconc  #ion concentrations
    profiles['nprodt']=nprodt #length
    profiles['kdiff_itemp']=kdiff_itemp  #variant flags
    profiles['kdiff_idens']=kdiff_idens
    profiles['variant']=variantid
    profiles['ion_temp']=ion_temp
    if rfxqlo_pro_variant:
        profiles['ion_dens']=ion_dens
        
    return profiles

#if profile_file:
#    toricnml['equidata']['idprof']= 1     
#    toricnml['equidata']['profnt_file']= profile_file     
#    toricnml['equidata']['nspec']= len(dt['iatm'])    
#    toricnml['equidata']['azi']  = dt['iazi']    
#    toricnml['equidata']['atm']  = dt['iatm']
#    toricnml['equidata']['aconc']= dt['aconc']

def toric_eqmodes(eq):
  #todo: get xmap and zmap from eq

  Xmap=eq.get('xmap') ; Zmap = eq.get('zmap')
  nmhd,ntheta=Xmap.shape ; imom=12
  rmc2d=np.zeros([nmhd,imom+1])
  rms2d=np.zeros([nmhd,imom+1])
  zmc2d=np.zeros([nmhd,imom+1])
  zms2d=np.zeros([nmhd,imom+1])

  for i in range(nmhd):
    cX = sft.fft(Xmap[i,:])/float(ntheta)
    cZ = sft.fft(Zmap[i,:])/float(ntheta)
    rmc2d[i,0] = np.real(cX[0])
    rmc2d[i,1:]= np.real(cX[1:imom+1]+np.flip(cX)[0:imom])
    rms2d[i,1:]= np.real((cX[1:imom+1]-np.flip(cX)[0:imom])*complex(0.,1,))
    zmc2d[i,0] = np.real(cZ[0])
    zmc2d[i,1:]= np.real(cZ[1:imom+1]+np.flip(cZ)[0:imom])
    zms2d[i,1:]= np.real((cZ[1:imom+1]-np.flip(cZ)[0:imom])*complex(0.,1,))

  eq['rzmcs2d'] = [rmc2d,rms2d,zmc2d,zms2d]


def write_equigs(eq,equigsfile):
    """ eg format
    Major radius (central)
    0.674548094E+00
    Major radius (axis)
    0.685074449E+00
    Magnetic field at major radius
    0.536006546E+01
    Total toroidal current
    0.685094066E+00
    Number of poloidal modes
    9
    Number of radial mesh points
    165
    Radial mesh

     write(ilun,'(A)')  'Fourier equilibrium coefficients'
     write(ilun,'(4E18.9)') (zrc(0,i),i=0,inx)
     write(ilun,'(4E18.9)') (zzc(0,i),i=0,inx)
     do  m=1,kmom
        write(ilun,'(4E18.9)') (zrc(m,i),i=0,inx)
        write(ilun,'(4E18.9)') (zzs(m,i),i=0,inx)
        write(ilun,'(4E18.9)') (zrs(m,i),i=0,inx)
        write(ilun,'(4E18.9)') (zzc(m,i),i=0,inx)
     enddo

    """
    equigs={}
    toric_eqmodes(eq)
    with open(equigsfile, "w") as file:
        torlheq_psimodes,torlheq_mmodes=eq['xmap'].shape
        psimap=eq['psipolmap']  #chosen uniform mesh

        file.write(' Major radius (central)(m) [generated by plasma.equigs.py]\n')
        file.write(f"{eq['rcentr']:18.9E}\n")
        equigs["rtorm"] =eq['rcentr']

        file.write(' Major radius (axis)(m)\n')
        file.write(f"{eq['rmaxis']:18.9E}\n")
        equigs["rmaxis"] =eq['rmaxis']

        file.write(' Magnetic field at major radius (m)\n')
        file.write(f"{eq['bcentr']:18.9E}\n")
        equigs["bcentr"] = np.abs(eq['bcentr'])
        equigs["sign_bcenter"] = np.sign(eq['bcentr'])

        file.write(' Total toroidal current (kA)\n')
        equigs["torcur"] =np.abs(eq['current']/1000.)
        equigs["sign_torcur"] =np.sign(eq['current']/1000.)
        file.write(f"{equigs['torcur']:18.9E}\n") #eqdsk is in Amps, torlh in kAmps

        equigs["rzmcs2d"]=eq['rzmcs2d']  #these are gotten from mapper.py and ffts
        rmc2d,rms2d,zmc2d,zms2d=eq['rzmcs2d']
        imom=rmc2d.shape[1]
        file.write(' Number of poloidal modes\n')
        equigs["imom"] = imom-1
        file.write(f"{equigs['imom']:5}\n")

        file.write(' Number of radial mesh points\n')
        equigs["nmhd"] = torlheq_psimodes
        file.write(f"{torlheq_psimodes:5}\n")

        file.write(' Radial mesh\n')
        rhopol=eq['rhopolmap']  #sqrt FluxGrid
        equigs["srad"] = rhopol
        formattedwrite(file,rhopol)

        file.write(' Fourier equilibrium coefficients\n')
        formattedwrite(file,rmc2d[:,0]) #dc modes
        formattedwrite(file,zmc2d[:,0])
        for i in range(1,imom):
            formattedwrite(file,rmc2d[:,i])
            formattedwrite(file,zms2d[:,i])
            formattedwrite(file,rms2d[:,i])
            formattedwrite(file,zmc2d[:,i])

        file.write(' Safety factor\n')
        qmap=np.interp(eq['psipolmap'],eq['fluxGrid'],eq['qpsi'])
        equigs["qqf"] = qmap
        formattedwrite(file,qmap)

        file.write(' Current profile [kA]\n')
        equigs["jcurr"]=np.abs(eq['Ipsi']/1000.0)
        formattedwrite(file,equigs["jcurr"])

        file.write(' Covariant B_phi, R*B_phi (m*T)\n')
        gmap=np.interp(eq['psipolmap'],eq['fluxGrid'],eq['fpol'])
        equigs["gcov"] = np.abs(gmap)
        formattedwrite(file,gmap)

        file.write(' Rho toroidal\n')
        equigs["rhotor"]=eq['rhotormap']
        formattedwrite(file,eq['rhotormap'])

        file.write(' Fraction Psi poloidal at last surface\n')
        equigs["lastpsi"]=1.0 #make better
        file.write(f"{equigs['lastpsi']:18.9E}\n")

    return equigs


def read_equigsfile(equigsfile='equigs.data'):
    "Read the equilibrium file created by toric in toricmode='equil',isol=0."

    def __get_varname(f, debug=False):
        "Reads next line from file f and returns it, optionally printing it."
        varname=f.readline()
        if debug:
            print (f.name,varname)
        return varname

    equigs_hdl=open(equigsfile,'r')
    equigs = {}
    equigs['file']=equigsfile
    
    varname = __get_varname(equigs_hdl)
    equigs["rtorm"] = np.fromfile(equigs_hdl,sep=" ",
                                  count=1,dtype=float)[0]

    varname = __get_varname(equigs_hdl)
    equigs["rmaxis"]= np.fromfile(equigs_hdl,sep=" ",
                                 count=1,dtype=float)[0]

    varname = __get_varname(equigs_hdl)
    equigs["bcentr"] = np.fromfile(equigs_hdl,sep=" ",
                                  count=1,dtype=float)[0]

    varname = __get_varname(equigs_hdl)
    equigs["torcur"]= np.fromfile(equigs_hdl,sep=" ",
                                  count=1,dtype=float)[0]

    varname = __get_varname(equigs_hdl)
    equigs["imom"] = np.fromfile(equigs_hdl,sep=" ",
                                 count=1,dtype=int)[0]
    imom = equigs["imom"]

    varname = __get_varname(equigs_hdl)
    equigs["nmhd"] = np.fromfile(equigs_hdl,sep=" ",
                                 count=1,dtype=int)[0]
    nmhd = equigs["nmhd"]

    varname = __get_varname(equigs_hdl)
    equigs["srad"] = np.fromfile(equigs_hdl,sep=" ",
                                 count=nmhd,dtype=float)

    #this needs to be reshaped or remapped into the R,Z sin cos
    #arrays toric uses
    varname = __get_varname(equigs_hdl)
    equigs["rzmcs2d"] = np.fromfile(equigs_hdl,sep=" ",
                                    count=2*nmhd+4*nmhd*imom,dtype=float)

    varname = __get_varname(equigs_hdl)
    equigs["qqf"] = np.fromfile(equigs_hdl,sep=" ",
                                count=nmhd,dtype=float)

    #logic checking for "END"
    varname = __get_varname(equigs_hdl)
    equigs["jcurr"] = np.fromfile(equigs_hdl,sep=" ",
                                  count=nmhd,dtype=float)
    
    varname = __get_varname(equigs_hdl)
    equigs["gcov"]= np.fromfile(equigs_hdl,sep=" ",
                                count=nmhd,dtype=float)

    varname = __get_varname(equigs_hdl)
    equigs["rhotor"] = np.fromfile(equigs_hdl,sep=" ",
                                   count=nmhd,dtype=float)

    varname = __get_varname(equigs_hdl)
    equigs["lastpsi"] = np.fromfile(equigs_hdl,sep=" ",
                                    count=1,dtype=float)[0]
   
    equigs_hdl.close()
    
    return equigs


def plot_equigs(equigs, ntheta=65, ax=None):
    if ax is None:
        ax = plt.gca()
        
    imom=equigs['imom'] 
    nmhd=equigs['nmhd']
    rzmcs2d=equigs['rzmcs2d']
    rmc2d0 =rzmcs2d[0:nmhd]  #center of each flux surface. first term is magnetic axis.
    zmc2d0 =rzmcs2d[nmhd:2*nmhd]
    rz=     rzmcs2d[2*nmhd:].reshape( (nmhd,4,imom), order='F' )

    Raxis=rmc2d0[0]
    rminor= (np.sum(rz[:,0,:],1)+ rmc2d0) - Raxis
    Rmajor=np.sum(rz[:,0,:],1) + rmc2d0

    #we keep the m=0 mode for sin coefficients for consistency
    rmc2d=np.zeros([nmhd,imom+1])
    rms2d=np.zeros([nmhd,imom+1])
    zmc2d=np.zeros([nmhd,imom+1])
    zms2d=np.zeros([nmhd,imom+1])

    rmc2d[:,0]=rmc2d0
    zmc2d[:,0]=zmc2d0
    rmc2d[:,1:]=rz[:,0,:]
    zms2d[:,1:]=rz[:,1,:]
    rms2d[:,1:]=rz[:,2,:]
    zmc2d[:,1:]=rz[:,3,:]

    rtest=np.zeros([nmhd,ntheta])
    ztest=np.zeros([nmhd,ntheta])
    theta = np.linspace(0,2.*np.pi,ntheta,endpoint=False)

    idx = np.linspace(0,imom+1,imom+1)
    for i in range(nmhd):
        for j in range(len(theta)):
            th=theta[j]
            rtest[i,j]=np.sum(rmc2d[i,:]*np.cos(th*idx)+ rms2d[i,:]*np.sin(th*idx) )
            ztest[i,j]=np.sum(zmc2d[i,:]*np.cos(th*idx)+ zms2d[i,:]*np.sin(th*idx) )

    #Plot coordinate mesh
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    #fig.set_figheight(12)

    #Theta lines
    for i in np.arange(0,len(rtest[0,:]),2):
        ax.plot(rtest[:,i],ztest[:,i])

    #Psi surface
    for i in np.arange(0,len(rtest[:,0]),5):
        ax.plot(rtest[i,:],ztest[i,:])
    ax.plot(rtest[-1,:],ztest[-1,:])

    ax.set_title('Toric Eq from equigs file '+equigs['file']);

    return ax


def XZ_from_equigs(equigs,dpsi=0,dtheta=0):
    "Not fully implemented, do not use"
    imom=equigs['imom'] 
    nmhd=equigs['nmhd']
    rmc2d0 =equigs['rzmcs2d'][0:nmhd]  #center of each flux surface. first term is magnetic axis.
    zmc2d0 =equigs['rzmcs2d'][nmhd:2*nmhd]
    rz=equigs['rzmcs2d'][2*nmhd:].reshape( (nmhd,4,imom), order='F' )
    rmc2d=rz[:,0,:]
    zms2d=rz[:,1,:]
    rms2d=rz[:,2,:]
    zmc2d=rz[:,3,:]
    return

    

def plot_diag(diag, eq0):
    #pgridx,pgridy,tgridx,tgridy,igsmhd,iqtest,ncopsi,jptheta,ntt,lpl=read_diag(diag)
    print ("toric.asc: Settings and dimensions: ",igsmhd,iqtest,ncopsi,jptheta,ntt,lpl)
    fig_grid=plt.figure()
    ax=fig_grid.add_subplot(111)
    ax.set_aspect('equal')
    fig_grid.suptitle('Eq from diag output',fontsize=16)
    for i in range(ncopsi):
        plt.plot(diag.pgridx[i],diag.pgridy[i],'k')
        for i in range(ntt):
            plt.plot(diag.tgridx[i],diag.tgridy[i],'b')
        #plt.plot(splotx,sploty,'g');
        plt.plot(diag.pgridx[0],diag.pgridy[0],'g');
        plt.contour(eq0['r']*100-42,eq0['z']*100,eq0['psizr'].T,20,colors='r');
    return
    
    
def read_diag(diagfile='toric.asc',idebug=False): #,linetype='color',ppp=1,savefig=False):
    "Parse the ASCII formatted diagnostic output for the equilibrium mesh. Not fully implemented."

    diag={}
    f0001=ff.FortranRecordReader('5i5')
    f0014=ff.FortranRecordReader('5i4')
    f0002=ff.FortranRecordReader('2i4')
    f0003=ff.FortranRecordReader('a10,1i4')
    f0004=ff.FortranRecordReader('5e16.9')
    f0005=ff.FortranRecordReader('a32')
    f0006=ff.FortranRecordReader('1e14.5')
    f0007=ff.FortranRecordReader('6e14.5')
    f0008=ff.FortranRecordReader('i5')
    f0009=ff.FortranRecordReader('6e13.4')
    

    with open(diagfile,'r') as of:
        [iqtest]=f0008.read(next(of))
        [plotName]=f0005.read(next(of))
        if idebug : print('plotName',plotName)
        [igsmhd, iudsym, npsi, ipsi, modmhd]=f0001.read(next(of))

        srad  = np.zeros(npsi, dtype=float)
        irad  = np.zeros(ipsi, dtype=float)
        rawcf = np.zeros(npsi, dtype=float)
        intcf = np.zeros(ipsi, dtype=float)

        
        def getpsiarray(of,fmt=f0009,nx=npsi):
            [plotName]=f0005.read(next(of))
            if idebug : print('plotName',plotName)
            return readArray(of,fmt,[nx],6)

        
        srad=getpsiarray(of,fmt=f0007)
        irad=getpsiarray(of,fmt=f0007)
#        [plotName]=f0005.read(next(of))
#        if idebug : print('plotName',plotName,npsi)
#        srad=readArray(f0007,[npsi])
#        [plotName]=f0005.read(next(of))
#        if idebug : print('plotName',plotName)
#        irad=readArray(f0007,[ipsi])

        if idebug : print('modmhd',modmhd)
        for m in range(modmhd+4):
            print('m',m)
            [plotName]=f0005.read(next(of))
            if idebug : print('plotName',plotName)
            rawcf=readArray(of,f0007,[npsi],6)
            intcf=readArray(of,f0007,[ipsi],6)

        qq=getpsiarray(of)
        aj=getpsiarray(of)
        ai=getpsiarray(of)
        vi=getpsiarray(of)
        ai2=getpsiarray(of)

        #Magnetic Field Configuration 
        [plotName]=f0005.read(next(of))                      
        if int(plotName)==0:
            print('idlout set to 0 in torica.inp. Stopping')
            return diag

        [ncy1, ncy2, nres, ncof, npvert]=f0014.read(next(of))
        aux_x=np.zeros(npvert, dtype=float)
        aux_z=np.zeros(npvert, dtype=float)

       # if ncy1:
            
        
def read_diag2(diagfile='toric.asc'):
    "Parse the ASCII formatted diagnostic output for the equilibrium mesh"
    import re #works but could be converted to fortran formatted read
    toric_diag=open(diagfile).readlines()

    gridstart=re.compile('Magnetic configuration')
    for idx in range(len(toric_diag)):
        l=toric_diag[idx]
        if gridstart.findall(l):
            lgrid=idx
            print (l=='Magnetic configuration\n',toric_diag[lgrid])

    #Get settings
    igsmhd,iqtest=np.fromstring(toric_diag[lgrid+1],dtype=int,sep=' ')
    ncopsi,jptheta=np.fromstring(toric_diag[lgrid+2],dtype=int,sep=' ')

    #formatting is 6E13.8, jptheta/6 should be number of lines
    nlines=int(np.ceil(jptheta/6))

    splotx=np.fromstring(''.join(toric_diag[lgrid+3:lgrid+3+nlines]),sep=' ')
    sploty=np.fromstring(''.join(toric_diag[lgrid+3+nlines:lgrid+3+2*nlines]),
                         sep=' ')

    idx=lgrid+3+2*nlines
    title=toric_diag[idx]
    print("title",title,nlines)
    idx+=1

    pgridx=[]
    pgridy=[]
    pgridx.append(splotx)
    pgridy.append(sploty)
    for i in range(ncopsi+1):
        pgridx.append(np.fromstring(''.join(toric_diag[idx:idx+nlines]),
                                    sep=' '))
        idx+=nlines
        pgridy.append(np.fromstring(''.join(toric_diag[idx:idx+nlines]),
                                    sep=' '))
        idx+=nlines

    title=toric_diag[idx]
    idx+=1

    ntt,lpl=np.fromstring(toric_diag[idx],dtype=int,sep=' ')
    idx+=1
    tgridx=[]
    tgridy=[]
    nlines=int(np.ceil(lpl/6))
    for i in range(ntt):
        tgridx.append(np.fromstring(''.join(toric_diag[idx:idx+nlines]),
                                    sep=' '))
        idx+=nlines
        tgridy.append(np.fromstring(''.join(toric_diag[idx:idx+nlines]),
                                    sep=' '))
        idx+=nlines

    #print ("Idx group is ",toric_diag[idx])
    diag={'pgridx':pgridx,'pgridy':pgridy,'tgridx':tgridx,'tgridy':tgridy,
          'igsmhd':igsmhd,'iqtest':iqtest,'ncopsi':ncopsi,'jptheta':jptheta,
          'ntt':ntt,'lpl':lpl}
    return diag #pgridx,pgridy,tgridx,tgridy,igsmhd,iqtest,ncopsi,jptheta,ntt,lpl




def stix_temperature(Prf,Te,ne,A,Z,Chi):
    """
    T effective from Stix 1975 for minority heating
    1.32e9*np.sqrt(3.14159)/(5.64e4**2*1.32e3**2)*2*np.sqrt(3.14159)*3.14159/20/9.11e-28* 1e7/1e28

    = 0.258 * (ln Lambda/20) . . .
    $$
    xi_{mathrm{mino}}^{mathrm{(Stix)}} approx{
    frac{0.258 (lnLambda/20) , [ T_{e}(mathrm{keV}) ]^{1/2}
             A_{mathrm{mino}} langle P_{mathrm{RF}} rangle_{mathrm{MW/m^{3}}}}
         {n_{e,20}^{2} , Z_{mathrm{mino}}^{2} , X_{mathrm{mino}}}
    }
    $$


    """

    lnlambda=24-np.log (np.sqrt(ne*1.E14)/ (Te*1000) ) #~21 for sparc Hmode

    xi =  0.258 * (lnlambda/20) * np.sqrt(Te)*A*Prf/( (ne *Z)**2 * Chi )

    return Te*(1+xi)


class toric_analysis:
    """
    Class to encapsulate tools used for toric analysis.
    depends on: python3, matplotlib, numpy, scipy, socket, time,
                         periodictable, f90ml, fortranformat
    
    Typical invocation:
    import toric_tools
    R=toric_tools.toric_analysis('toric.ncdf',mode='ICRF') #'LH' for lower hybrid
          toric_data keyword typically one of: 'toric.data', 'toric_cfg.nc', fort.9'
    R.plot_2Dfield(component='Re2Ezeta',logl=10) #etc

    Important functions:
    R.info() # netcdf contents and metadata
    R.plotpower(power='PwIF',species=1) #different power profiles
    R.plot_2Dfield(component='Re2Ezeta',logl=10) #Two Dim plots of quantities
    R.plot_1Dfield(component='Re2Ezeta') #One Dim plots of quantities
    R.threeplots() #2D Ez, electron power and poynting flux and poloidal
                   #spectrum. Now additional 2D power by species so more than three
    """


    def __init__ (self, toric_name='toric.ncdf', toric_data="toric.data",
                  mode='ICRF', idebug=False, comment='', layout='poster',
                  path="./"):
        import socket
        from time import gmtime

        self.toric_name=toric_name
        self.toric_data=toric_data
        self.mode = mode
        self.__version__ = 1.1
        self.idebug = idebug

        self.mylw=1.0
        self.mypt=18.0
        self.fsc=4.0
        self.fw='bold'
        self.set_layout(layout)

        self.path=path

        self.prov = {"user":"noname","host":"noname","gmtime":"notime","runid":"noid",
                "path":"", "comment":""}
        self.label = True
        self.equigs = {}
        self.toricdict={}
        self.nml=f90nml.read(os.path.join(path,'torica.inp') )

        if (self.mode[:2]=='LH'):
            self.namemap={'xpsi':'tpsi','poynt':'vpoynt','pelec':'S_eld',
                 'e2d_z':'E2d_z_re','xplasma':'x_plasma', 'zplasma':'z_plasma',
                          'xeqpl':'xeqpl'}
            if self.toric_name=='None': self.toric_name='TORICLH.cdf'
        else:
            self.namemap={'xpsi':'Pw_abscissa','poynt':'PoyFlx','pelec':'PwE',
                         'e2d_z':'Re2Ezeta','xplasma':'Xplasma',
                          'zplasma':'Zplasma', 'xeqpl':'Ef_abscissa'}

##Open the toric netcdf file
        try:
            self.cdf_hdl = netcdf_file(path+self.toric_name,mmap=False )
            dvs = self.cdf_hdl.variables
        except IOError:
            print ('CRITICAL: ',self.toric_name,' not found.')
            self.cdf_hdl = None
            return 

        try:
            self.qlde_hdl = netcdf_file(path+"toric_qlde.cdf",mmap=False)
        except IOError:
            print ('Non-CRITICAL: ',path+"toric_qlde.cdf",' not found.')
            self.qlde_hdl = None

        try:
            self.data_hdl = netcdf_file(path+self.toric_data,mmap=False )
        except IOError:
            print ('CRITICAL: ',self.toric_data,' not found.')
            self.data_hdl = None

        xx = dvs[self.namemap['xplasma']].data
        nant=1
        self.nspec=self.nml['equidata'].get('nspec')
        if not self.nspec: print('Warning nspec not found in namelist')

        self.spec=self.get_spec()
        if self.data_hdl:
            dv=self.data_hdl.variables
            self.nspec=self.data_hdl.dimensions['n_of_ionspec']
            ant_ipsi= (np.abs(xx[0,:] - dv['antenna_radius'].data[0])).argmin()
            
            self.antenna={
                'nant':nant, 'length':dv['ant_length'].data[0],
                'theta':dv['ant_position'].data[:nant],
                'rmajor':dv['antenna_radius'].data[0]+dv['axis_radius'].data[0],
                'radius':dv['antenna_radius'].data[0],'ipsi':ant_ipsi }
        else:
            self.antenna={
                'nant':1, 'length':10,'theta':0.,
                'rmajor':dvs['Raxis'].data+dvs['xedg_out'].data,
                'radius':dvs['xedg_out'].data,'ipsi':1 }

        self.prov["host"]=socket.getfqdn()
        self.prov["user"]=os.getenv("USER")
        self.prov["gmtime"]=gmtime()
        self.prov["comment"]=comment
        self.prov["path"]= os.path.abspath('')+'/'+self.toric_name

#marker sequence
        self.markers = ['o','1','2','s','*','+','H','x']
        self.ls = ['-','--','-.',':','--','-.',':']
        self.bw = False

        return


    def close (self):
        try:
            self.cdf_hdl.close()
        except IOError:
            print ('CRITICAL: ',self.toric_name,' not found.')

        try:
            self.qlde_hdl.close()
        except IOError:
            print ('Non-CRITICAL: ',path+"toric_qlde.cdf",' not found.')

        try:
            self.data_hdl.close()
        except IOError:
            print ('Non-CRITICAL: ',path+self.toric_data,' not found.')

        return


    def info( self ):
        "Prints a list of the contents of present TORIC3D output files"

        if self.cdf_hdl:
            for hdl in [self.cdf_hdl]:
                print ('The toric file, ',self.toric_name,', contains:')
                print ('----------------------------------------------')
                print ("The global attributes: ",self.cdf_hdl.dimensions.keys())
                print ("File contains the variables: ", self.cdf_hdl.variables.keys())

        if self.qlde_hdl:
            for hdl in [self.qlde_hdl]:
                print ('The toric file, ',"toric_qlde.cdf",', contains:')
                print ('----------------------------------------------')
                print ("The global attributes: ",self.qlde_hdl.dimensions.keys()  )
                print ("File contains the variables: ", self.qlde_hdl.variables.keys())

        if self.data_hdl:
            for hdl in [self.data_hdl]:
                print ('The toric file, ',self.toric_data,', contains:')
                print ('----------------------------------------------')
                print ("The global attributes: ",self.data_hdl.dimensions.keys()  )
                print ("File contains the variables: ", self.data_hdl.variables.keys())
              #  print ("antenna", self.data_hdl.var

        print ('----------------------------------------------')
        print ("Provenance metadata: ", self.prov)

        #print ('----------------------------------------------')
        #print ('TORIC parameters')

        print ('----------------------------------------------')
        print ('Power partitions')
        print ('Species                 '+
               ''.join([spec.get('name')+'  ' for spec in self.get_spec()]) )
        for var in ['TPwIF', 'TPwIH', 'TPwEFW', 'TPwEIBW']:
            print("{0:} {1:>3} ".format(self.cdf_hdl.variables[var].
                    long_name.decode('UTF-8') ,
                    ListToFormattedString(self.cdf_hdl.variables[var].data,'{:.2f}%') ) )
        print ('----------------------------------------------')
        return


    def get_spec( self ):
        "Collect info on species for ICRF sim in TORIC in nice readable format"

        spec_toric=[]
        if self.data_hdl:
            dv=self.data_hdl.variables
            atm = dv['atomic_masses'][:]
            atz = dv['atomic_charges'][:]
            conc= dv['concentrations'][:]
            for i in range(len(atm)):
                name=str(elements[atz[i]][atm[i]])
                spec_toric.append({ 'name':name,'A':atm[i],'Z':atz[i],
                               'Conc%':100*conc[i] })
            spec_toric.insert(0,{'name':'e', 'A':0, 'Z':-1 , 'Conc%': 100})
        else:
            print('Warning no species info found')
            for i in range(8):            
                spec_toric.append({'name':'none','A':1,'Z':1,'Conc%':100})
            spec_toric.insert(0,{'name':'e', 'A':0, 'Z':-1 , 'Conc%': 100})
        return spec_toric

                  
    def plotb0( self, ir=45, db=0, eps=0 ):
        """Contour plot of bounce averaged Dql coefficient, dB0 """
        if not self.qlde_hdl:
            print ("qlde File not found")
            return
        if (eps>0):
            db=2*eps/(1.-eps)

        dqlpsi=self.qlde_hdl.variables['Psi'].data
        dqltemp=self.qlde_hdl.variables['Tem'].data
        dql_LD=self.qlde_hdl.variables['Qldce_LD'].data
        nuperp=self.qlde_hdl.dimensions['VelPrpDim']
        nupar=self.qlde_hdl.dimensions['VelDim']

        umax=(self.qlde_hdl.variables['Umax'].data)[0]
        umin=(self.qlde_hdl.variables['Umin'].data)[0]
        upar=np.arange(nupar)/float(nupar-1)*(umax-umin)+umin
        uperp=np.arange(nuperp)/float(nuperp-1)*umax
        vx,vz=np.meshgrid(uperp,upar)

        fig=plt.figure(figsize=(2.*8.3,2.*3.7))
        plt.axes().set_aspect(1, 'box')
        #plot passing trapped boundary
        roa = dqlpsi[ir]
        if (eps<0):
            db=2.*np.abs(eps)*roa/(1.-np.abs(eps)*roa)

        if (db>0):
            vpar=np.sqrt(db)*umax
            plt.plot([0,vpar],[0,umax],'k',[0,-vpar],[0,umax],'k',linewidth=2)

        dq=np.transpose(np.log((dql_LD[ir,:,:])+1.)/np.log(10)) #np.abs
        mxdq=int(dq.max())
        ll=range(mxdq-10,mxdq)
        cd=plt.contourf(vz,vx,dq,levels=ll)
#,10)
        plt.gca().set_ylim(0,umax)
        cbar=plt.colorbar(cd)

        plt.title(r'log10 $\lambda$<B> at r/a='+str(roa)[0:4],size=30)
        plt.ylabel(r'$u_{\bot0}/u_{n}$',size=20)
        plt.xlabel(r'$u_{||0}/u_{n}$',size=20)
        plt.draw() #make sure ylimits are updated in displayed plot
        plt.close()
        return


    def plotpower( self, xaxis=None, power=None, species=None ):
        """
        Plot power profiles versus specified radius for all, or
        listed species. Overplots by default. species is 0 indexed.
        """

        l=-1
        if (self.mode[:2]=='LH'):
            if (xaxis==None):
                xaxis='tpsi'
            l=self.__plot1D(xaxis,'S_eld','Power absorbed on electrons')
            plt.xlabel(r'$\sqrt{\psi_{pol}}$')
        else:
            if (xaxis==None):
                xaxis='Pw_abscissa'
            if (power==None):
                power='PwE'
            if (species==None): # or species==0):
                l=self.__plot1D(xaxis,power)
            else:
                nspec=self.cdf_hdl.dimensions['SpecDim']
                if (species<=nspec and species>0):
                    l=self.__plot1D(xaxis,power,idx2=species-1)#zero indexing
                else:
                    print("Invalid species label:"+str(species))
            plt.xlabel('X/a')
        cf=plt.gcf()
        cf.subplots_adjust(bottom=0.14)

        return l


    def psiplot( self, y ):
        """
        "Plot versus rhopsi. Returns handle on line to modify line
        style if desired using setp.
        """
        psi=self.namemap['xpsi']

        line=self.__plot1D(psi,y)
        plt.xlabel(r'$\sqrt{\psi_{pol}}$')
        plt.ylabel(y)

        return line


    def plot_1Dfield( self, component ):
        "Field versus midplane specified."

        line=self.__plot1D(self.namemap['xeqpl'],component,
                           'Wave Field component on the midplane',component)
        plt.xlabel(r'$X[cm]$')
        return line


    def __plot1D( self, xvar, yvar, ptitle=None, plabel='', idx2=None):
        "Internal 1D plot"
        x=self.cdf_hdl.variables[xvar]
        y=self.cdf_hdl.variables[yvar]
        if (self.mode[:2]=='LH'):
            xname=''
            yname=plabel
        else:
            xname=x.long_name.decode('UTF-8')
            yname=(y.long_name[0:20]).decode('UTF-8') + y.units.decode('UTF-8')

        x=x[:]
        if (idx2!=None):
            if (len(y.shape)==2):
                y=y[:,idx2]
            else:
                print(yvar+" has wrong number of dims in __plot1d.")
                return
        else:
            y=y[:]

        if (np.size(y) > np.size(x)):
            print ("ToricTools.__plot1D resizing",yvar)
            y=np.array(y)[0:np.size(x)]

        line=plt.plot(x,y)
        if ptitle:
            plt.title(ptitle)
        plt.xlabel(xname)
        plt.ylabel(yname)


        return line


    def __map_celef (self):
        """
        Mapping toric field in celef.cdf to re and im parts of the three
        components for LH mode.
        """

        return


    def __getvar__( self, name ):
        """
        Internal function to retrieve variable from data file with checking.
        """

        try:
            value=self.cdf_hdl.variables[name].data
        except NameError:
            print ('CRITICAL: variable not found')
            exit #raise Exception,'Bad variable name in getvar: %s' % name

        return value


    def fft( self, component='undef',maxr=1. ):

        if (self.mode[:2]=='LH'):
            radius = 'psime'
        if (component=='undef'):
            if (self.mode[:2]=='LH'):
                component='E2d_z_re'
            else:
                component='Re2Ezeta'

        field = self.__getvar__(component)
        rad   = self.__getvar__(radius)
        #field taken to be 2D with shape (ntheta,npsi)
        ntt=field.shape[0]
        nelm=int(field.shape[1]*maxr)
        nlevels=100
        levels=np.arange(nelm/nlevels,nelm-1,nelm/nlevels)
        fftfield = np.zeros((ntt,levels.shape[0]),'complex128')
        i=0
        for ir in levels:
            ffield = (ft.fft(field[:,ir]))
            fftfield[:,i] = ffield
            i=i+1

        return fftfield


    def spectrum( self, component='undef',maxr=1.,cx=0,levels=-1, q=None,
                  figname=None):
        """Calculate poloidal spectrum of two dimensional field component.
        """

        if (self.mode[:2]=='LH'):
            radius = 'psime'
        else:
            radius = 'Pw_abscissa'

        if (component=='undef'):
            if (self.mode[:2]=='LH'):
                component='E2d_z_re'
                componenti='E2d_z_im'
            else:
                component='Re2Eplus'
                componenti='Im2Eplus'

        f=plt.figure()

        if (component=="power"):
            field = self.get_power2D()
        else:
            field = (self.__getvar__(component))#[:,:]

        if (cx==1):
            fieldi = (self.__getvar__(componenti))#[:,:]
            field=np.array(field)+1.0j*np.array(fieldi)

        rad   = self.__getvar__(radius)

        #field taken to be 2D with shape (ntheta,npsi)
        field=field+1.e-20
        ntt=field.shape[0]
        #nelm=int(field.shape[1]*maxr)
        nelm=int(np.size(rad)*maxr)
        if self.idebug : print(levels)
        if (np.size(levels)==1):
            nlevels=7
            levels=(np.arange(nlevels)*nelm*1./nlevels).astype(int)
        else:
            levels=(np.array(levels)*nelm).astype(int)
            nlevels=np.size(levels)

        levels=levels[1:nlevels]
        rlevels=rad[levels]
        if self.idebug: print('SPECTRUM of ', component, nelm,nlevels,levels,rlevels)

        th = np.arange(ntt)-ntt/2

        ymax = 0.
        ymin = 0.

        i=0
        thq=th
        for indr in range(levels.size): #levels:
            #fft in python isn't normalized to N
            ir=levels[indr]

            if q!=None: #JCW fix
                thq=-2.5*(1+0.3)/(1+0.3*rlevels[indr])*(1+th/191./q(rlevels[indr]))

            ffield = ft.fftshift(np.log10(abs(
                ft.fft(field[:,ir]))/float(ntt)+1.e-20))
            ymax = np.max( [ymax, np.max(ffield)] )
            ymin = np.min( [ymin, np.min(ffield)] )
            plabel='%5.2f' % rad[ir]
            if self.bw:
                plt.plot( thq, ffield, label=plabel,
                          linestyle=self.ls[i],color='k')
                i=i+1
            else:
                plt.plot( thq, ffield, label=plabel )

        ffield = ft.fftshift(np.log10(abs(
            ft.fft(field[:,nelm-1]))/float(ntt)+1.e-5))
        ymax = np.max( [ymax, np.max(ffield)] )
        ymin = np.min( [ymin, np.min(ffield)] )

        #plot antenna spectrum
        plabel='ant'
        if self.idebug: print ("range, levels", rlevels)
        if self.idebug: print ("ymax", ymax,ymin)
        plt.plot( thq, ffield,  label=plabel, color='grey')

        if q!=None:
            plt.axis( xmin=-8,xmax=8 )
        else:
            plt.axis( xmin=-ntt/4, xmax=ntt/4)

        plt.axis( ymin=-10)
        plt.legend(loc=(1.05,0),title='r/a surface',labelspacing=.1,fontsize=10)
        plt.xlabel('m')
        plt.ylabel('log10 scale')
        plt.title('Poloidal spectrum')
        plt.ylim(bottom=-5)
        plt.tight_layout()
        if figname:
            plt.savefig(figname+'.pdf',format='pdf')
            plt.savefig(figname+'.png',format='png')
        plt.close()
        return


    def set_layout( self, layout='poster' ):

        if (layout == 'paper'):
            self.mylw=2.0
            self.mypt=10.0
            self.fsc=2.0
            self.fw='normal'

        if (layout == 'poster'):
            self.mylw=3.0
            self.mypt=20.0
            self.fsc=4.0
            self.fw='bold'


        params = {
            'axes.linewidth': self.mylw,
            'lines.linewidth': self.mylw,
            'axes.labelsize': self.mypt,
            'font.size': self.mypt,
            'legend.fontsize': self.mypt,
            'axes.titlesize': self.mypt+4.0,
            'xtick.labelsize':self.mypt-2,
            'ytick.labelsize':self.mypt-2,
            'font.weight'  : self.fw,
            'text.usetex' : False
            }
        plt.rcParams.update(params)

        return


    #note that if plot commands are in the toplevel, they will not return
    #to the prompt, but wait to be killed.
    def plot_2Dfield(self, component='E2d_z',species=None,logl=0,
                     xunits=1.0,axis=(0.0,0.0), im=False, cmap=None,
                     scaletop=1.0,scalebot=1.0,ax='undef',fig='undef',
                     maxsurface=0.99,lscaletop=0.0,lscalebot=0.0):
        """
        Example of using netcdf python modules to plot toric solutions
        requires numpy and matplotlib and netcdf modules for python.

        To overplot with limiter, made from efit plotter:
        R.plot_2Dfield(component='Im2Eplus',logl=20,xunits=0.01,
                       axis=maxis,fig=fig1)

        Easier is to plot solution first, then overplot limiter,
        scaled appropriately:
        p.plot ( rlim*100.-maxis[0], zlim*100.-maxis[1], 'k', linewidth = 2 )

        """
        if self.idebug: print('call args',locals() )
        #Select color table by log, Power, or field
        CT='jet'

        R0=axis[0]
        Z0=axis[1]
        barfmt='%5.2e' #'%4.1e' #'%3.1f'
        #what should colorbar with be? format=4.1e means 8 characters
        #the bar and title of the bar add about 4 characters.
        #there are 72.27 pt/in
        #12 characters * self.mypt /72.27 pt/in = #in
        legend_frac=12*self.mypt/72.27
        title=component

        xx  = self.cdf_hdl.variables[self.namemap['xplasma']].data
        yy  = self.cdf_hdl.variables[self.namemap['zplasma']].data

        if self.idebug: print (self.mode,'mode')
        if (self.mode[:2]=='LH'):
            if (im):
                im_e2dname=component+'_im'
                title='|'+component+'|'
                component=component+'_re'
        else:
            if species:
                title=title+' for '+ self.spec[species].get('name')
            if (component=='E2d_z'):
                component='Ezeta'

            if component[0]!='T':
                im_e2dname='Im2'+component
                if (im) : title='|'+component+'|'
                component='Re2'+component


        if (component=="power" and self.mode[:2]=='LH'):
            e2d = self.get_power2D()
        else:
            e2d = (self.cdf_hdl.variables[component]).data

        if (im):
            im_e2d=(self.cdf_hdl.variables[im_e2dname]).data
            e2d = abs(e2d+1.0j*im_e2d)

        if (self.mode[:2]!='LH' and species):
            if self.idebug: print('plot2D, indexing species', species)
            e2d = e2d[:,:,species-1]


        if self.idebug: print ("2D Matrix shape:", np.shape(xx))


    #contour with 3 args is confused unless arrays are indexed slices
    #need wrapper to close periodicity in theta direction for this
    #tricky, array indexing different from ncdf slicing
    #this step is needed because periodic dimension is not closed.
    #i.e. its [0,pi) not [0,pi]
        dd=np.shape(xx)
        sx=dd[0] #theta
        sy=dd[1] #psi
        lastpsi=int(sy*maxsurface)
        if (self.idebug): print("2D plot shapes:",sx,sy,lastpsi,maxsurface)

        xxx=np.zeros((sx+1,sy),'d')
        xxx[0:sx,:]=xx[:,:]
        xxx[sx,:]=xx[0,:]
        yyy=np.zeros((sx+1,sy),'d')
        yyy[0:sx,:]=yy[:,:]
        yyy[sx,:]=yy[0,:]

        xxx=(xxx+R0)*xunits
        yyy=(yyy+Z0)*xunits

        ee2d=np.zeros((sx+1,sy),'d')
        ee2d[0:sx,:]=e2d[:,:]
        ee2d[sx,:]=e2d[0,:]

        emax=np.max(ee2d[:,:lastpsi].ravel())
        emin=np.min(ee2d[:,:lastpsi].ravel())

        #contouring levels
        rmax=max([abs(emax),abs(emin)])*scaletop
        rmin=min([0.,emax,emin])*scalebot
        #val=arange(emin,emax,(emax-emin)/25.,'d')
        if self.idebug: print("2D rmax", rmax)
        if not rmax: rmax=1e4
        val=np.arange(-rmax*1.1,rmax*1.1,(rmax+rmax)/50.,'d')
        if (im):
            val=np.arange(rmin,rmax*1.1,(rmax)/24.,'d')
        if self.idebug: print ("values",val)

        #reverse redblue map so red is positive
           # revRBmap=cmap_xmap(lambda x: 1.-x, cm.get_cmap('RdBu'))

           #finally, make the plot
        cwidth=xxx.max()-xxx.min()
        cheight=yyy.max()-yyy.min()
        asp=cheight/cwidth
        if self.idebug: print ("plot aspect ratio:", asp)

        #leave space for bar
        if (fig=='undef'):
            fig=plt.figure(figsize=(self.fsc*3.0+legend_frac,3.0*self.fsc*asp))
            fig.subplots_adjust(left=0.02,bottom=0.15,top=0.90)

        sax=plt.axes().set_aspect(1, 'box')

        maxpsi=xxx.shape[1]-1
        plt.plot(xxx[:,maxpsi],yyy[:,maxpsi],'k-')

        #add LCF
        lcfpsi=self.cdf_hdl.dimensions['PsiPwdDim']
        plt.plot(xxx[:,lcfpsi],yyy[:,lcfpsi],'grey')


        #read ant length.  Calculate arc length vs theta to this value/2
        #in each direction, this plots the antenna location
        anthw=max(int(sx*0.01),4)
        ant_it_height= int(sx*self.antenna['length']/2/
                           ( 2.*np.pi * self.antenna['radius'] ) )

        ant_it_pos   =int(self.antenna['theta']*sx/360.)
        r1=np.arange(  ant_it_pos, ant_it_pos+ant_it_height+1)%sx
        r2=np.arange( (ant_it_pos-ant_it_height), (ant_it_pos+1))%sx
        plt.plot(  xxx[ r1, self.antenna['ipsi'] ],
                   yyy[ r1, self.antenna['ipsi'] ],
                   'orange',linewidth=4 )
        plt.plot(  xxx[ r2, self.antenna['ipsi'] ],
                   yyy[ r2, self.antenna['ipsi'] ],
                   'orange',linewidth=4 )

        if self.idebug: print("antenna: ", yyy[sx-anthw+1:sx+1,maxpsi], 'it:',
                              ant_it_pos,ant_it_height,sx)
        if self.label:
            ax=plt.gca()
            sublabel=self.prov['path']
            if self.idebug: print ('sublabel: ',sublabel)
            plt.text(-0.2,-0.3,sublabel,transform = ax.transAxes,fontsize=4)

        if (logl > 0):
            title='log10 '+title
            barfmt='%6.2e'

            ##labels and titles
            #xlabel(getattr(xx,'long_name')+'('+getattr(xx,'units')+')')
            #ylabel(getattr(yy,'long_name')+'('+getattr(yy,'units')+')')
            #title(getattr(e2d,'long_name')+'('+getattr(e2d,'units')+')')
        plt.xlabel('X(cm)')
        plt.ylabel('Z(cm)')
        plt.title(title,fontsize=self.mypt+2.0)


        if (logl <= 0):
            if not cmap: cmap='jet'
            CS=plt.contourf(xxx[:,:lastpsi],yyy[:,:lastpsi],
                            ee2d[:,:lastpsi],val,cmap=cmap) #cm.jet)
            #            for it in range(sx):
            #                PS=plt.plot(xxx[it,:lastpsi],yyy[it,:lastpsi],'k')
            #            for ip in range(sy):
            #                PS=plt.plot(xxx[:,ip],yyy[:,ip],'b')


        if (logl > 0):
            if not cmap: cmap='hot'
            lee2d=np.log(np.abs(ee2d[:,:lastpsi])+1.0)/np.log(10)
            rmax=lee2d.ravel()[lee2d[:,:lastpsi].argmax()]+lscaletop
            rmin=lee2d.ravel()[lee2d[:,:lastpsi].argmin()]+lscalebot

            if rmin==rmax:
                print('Warning rmin=rmax, no power,skipping',rmin,rmax,component,species)
                return
            val=np.arange(rmin,rmax,(rmax-rmin)/(logl*1.0),'d')
            CS=plt.contourf(xxx[:,:lastpsi],yyy[:,:lastpsi],
                            lee2d[:,:lastpsi],val,cmap=cmap) #cm.jet)

            ##put the contour scales on the plot
            #tricky, fraction needs to be specified to be part by which
            #horizontal exceed vertical

        cbar=plt.colorbar(CS,format=barfmt,fraction=0.05,pad=0.02) #ax=sax)
        cbar.ax.set_ylabel('levels')
        plt.tight_layout()


        if self.idebug: print ("contour values",CS.levels,'|',rmax,rmin)

        return CS,cbar


    ### user routines using the above, could be in a different module
    def get_power1D( self, species ):
        "Plots power across the midplant by averaging over Z"
        
        from scipy.interpolate import griddata
        
        XX=self.cdf_hdl.variables['Xplasma'][:]
        ZZ=self.cdf_hdl.variables['Zplasma'][:]
        X1D=self.cdf_hdl.variables['Ef_abscissa'][:]
        Z1D=np.linspace(np.min(ZZ),np.max(ZZ), len(X1D))
        XXcart,YYcart    = np.meshgrid(X1D,Z1D)

        pwr=self.cdf_hdl.variables['TDPwE']

        grid_e = np.reshape(griddata( ( XX[:,:].ravel(),ZZ[:,:].ravel() ),
                             pwr[:,:].ravel(),
                           ( XXcart.ravel(), YYcart.ravel() ),
                             method='nearest'
                                     ), (len(X1D),len(Z1D)) )
        slabpwr = np.sum(grid_e[:,:],axis=0)/float(len(Z1D))
        pwr1d = griddata( ( XX[:,:].ravel(),ZZ[:,:].ravel() ),
                             pwr[:,:].ravel(),
                           ( X1D[None,:],(X1D*0.)[None,:] ),
                             method='nearest'
                          ).T
        
        return X1D, pwr1d #grid_e, slabpwr, pwr1d

    
    def powpoynt( self ):
        "Plots powers and poynting flux"
        import matplotlib.colors as mcolors
        pcolors=list(mcolors.BASE_COLORS)

        fig = plt.figure(figsize=(12,9) )
        ax1 = fig.add_subplot(111)
        ax1.set_prop_cycle(color=['red', 'purple','orange', 'black',
                                  'green', 'blue', 'grey', 'gold', 'darkgreen'],
                  marker=['o', '+', 'x', 'o', 'v', '^', '<', '>', '.'])

        line1a,=self.psiplot(self.namemap['pelec'])
        #can use setp(lines, ) to change plot properties.
        plt.setp(line1a, label='electrons')
        line1b,=self.psiplot('PwEIBW')
        plt.setp(line1b, label=r'$P_{ibw}$')
        lines=[line1a,line1b]

        #add first two species if ICRF, add logic to plot if power percent
        #is larger than 0.5%
        if (self.mode[:2]!='LH'):
            nspec=self.cdf_hdl.dimensions['SpecDim']
            if self.idebug: print(self.nml['equidata'] )
            spec=self.get_spec()
            tpowerF=self.cdf_hdl.variables['TPwIF']
            tpowerH=self.cdf_hdl.variables['TPwIH']

            for ispec in range(nspec):
                if self.idebug:
                    print('total powers',tpowerF[ispec],tpowerH[ispec],ispec)
                if tpowerF[ispec]>0.1:
                    ltemp,=self.plotpower(power='PwIF',species=ispec+1)
                    plt.setp(ltemp,label='Fund '+spec[ispec+1]['name'])
                    print('Fund '+spec[ispec+1]['name'])
                    lines.append(ltemp)

                if tpowerH[ispec]>0.1:
                    ltemp,=self.plotpower(power='PwIH',species=ispec+1)
                    plt.setp(ltemp,label='Harm '+spec[ispec+1]['name'])
                    print('Fund '+spec[ispec+1]['name'])
                    lines.append(ltemp)



        ax2 = ax1.twinx()
        line20,=self.psiplot(self.namemap['poynt'])
        #line20.set_color('g')
        plt.setp(line20,label='<ExB>')

        lines.append(line20)
        if (self.idebug):
            for ll in lines:
                print('lines',ll, type(ll) ) #.get_label() )

        #set axis floor at 0
        ymax=np.average(self.cdf_hdl.variables[ self.namemap['pelec'] ][14:])*20
        #pnt('ymax',ymax,self.cdf_hdl.variables[ self.namemap['pelec'] ][:])

        ax1.set_ylim(0,0.5)
        ax2.set_ylim(0)
        ax1.set_ylabel('Power',color='b')
        
        #change color and symbol
        plt.setp(line20,color='r', label='<ExB>')
        ax2.set_ylabel('Poynting',color='r')
        ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

        #make  legend too
        plt.legend( handles=lines, loc='center right', 
                     ncol=1, fancybox=True, shadow=True)

        plt.xlim( 0, 1 )
        plt.tight_layout()
        plt.draw()
        return


    def powerion( self ):
        "Plots electron power and poynting flux"
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        line1,=self.psiplot(self.namemap['pelec'])
        #can use setp(lines, ) to change plot properties.
        plt.setp(line1,color='b',marker='+',label='seld')
        ax1.set_ylabel('Power_e',color='b')

        ax2 = ax1.twinx()
        line2,=self.psiplot(self.namemap['poynt'])
        #set axis floor at 0
        plt.gca().set_ylim(0)
        #change color and symbol
        plt.setp(line2,color='r',marker='.',label='Poynt')
        ax2.set_ylabel('Poynting',color='r')
        #make  legend too
        plt.legend( (line1,line2), (r'$P_{eld}$','<ExB>'),loc=2 )
        plt.axes().set_aspect(1, 'box')
        fig.subplots_adjust(left=0.12,bottom=0.12,top=0.96,
                            right=0.82,hspace=0.32)
        sax=plt.axes().set_aspect(1, 'box')

        plt.draw()
        return fig


    def xpsi_map( self ):
        """Return map of x(theta=0)/x(psi=1,theta=0) versus psipol."""
        xmap=1.0

        return xmap


    def get_power2D( self ):
        #figure out a sed way of cutting these lines into the file.
        #also need to replace '-0.' with ' -0.'
        #sed -n -e '/elec/,/,/p' filename | sed -e '/-0\./ -0./g' > torica_2dpower.sol
        try:
            toricsol = open('torica_2dpower.sol','r')
        except IOError:
            print ('CRITICAL: torica_2dpower.sol not found.')
            print ('Try to generate:')

            if (self.mode[:2]=='LH'):
                cmd="sed -n  '/elec/,$p' torica.sol| sed 's/-0\\./ -0./g' > torica_2dpower.sol"
            else:
                cmd="sed -n  '/elec/,$p' toric.sol| sed 's/-0\\./ -0./g' > torica_2dpower.sol"

            os.system(cmd)
            toricsol = open('torica_2dpower.sol','r')

        #skip title and max value
        toricsol.readline()
        toricsol.readline()

        if (self.mode[:2]=='LH'):
            nt=self.cdf_hdl.dimensions['ntt']
            nr=self.cdf_hdl.dimensions['mptpsi']
        else:
            nt=self.cdf_hdl.dimensions['ThetaDim']
            nr=self.cdf_hdl.dimensions['PsiPwdDim']

        power=np.fromfile(toricsol,sep=" ",count=nt*nr,dtype=float)
        toricsol.close()

        power=np.transpose(np.reshape(power,(nr,nt)))
        return power


    def threeplots( self, prefix='' ):
        """
        Makes and saves the three most commonly used plots. Plots are saved in
        the current directory. An optional prefix can be used to label them or
        change the save path.
        * Power and poynting flux on one plot as eps.
        * The polodial power spectrum on six flux surfaces for convergence
          as eps.
        * And the 2D parallel electric field contour plot as a png.
        """

        f1=plt.figure()

        self.spectrum(cx=1,levels=np.linspace(0.,0.98,12),figname=prefix+'spectrum')

        if (self.mode[:2]!='LH'):
            f2a=plt.figure()#figsize=(8,12))
            self.plot_2Dfield(component='Eplus', maxsurface=0.93,lscalebot=1,
                              lscaletop=0,im=True, logl=25,fig=f2a)
            plt.draw()
            plt.savefig('log10Eplus2d.png',format='png')

            f2b=plt.figure(figsize=(8,12))
            self.plot_2Dfield(component='Eplus', maxsurface=0.93,fig=f2b)
            #,scaletop=.4,scalebot=0.2)
            plt.draw()
            plt.savefig('Eplus2d.png',format='png')
            plt.close()

            f2b=plt.figure(figsize=(8,12))
            self.plot_2Dfield(component='Eminus', maxsurface=0.93,fig=f2b)
            #,scaletop=.4,scalebot=0.2)
            plt.draw()
            plt.savefig('Eminus2d.png',format='png')
            plt.close()            

        f3=plt.figure()#figsize=(8,12))
        self.plot_2Dfield(im=True,logl=25,fig=f3)#,scaletop=0.8) #default to Ez
        plt.draw()
        plt.savefig(prefix+'log10Ez2d.png',format='png')
        plt.close()
        
        f3=plt.figure()
        self.powpoynt()
        plt.draw()
        plt.savefig(prefix+'powerpoynt.pdf',format='pdf')
        plt.savefig(prefix+'powerpoynt.png',format='png')
        plt.close()
        
        f3=plt.figure()#figsize=(8,12))
        self.plot_2Dfield(component='TDPwE',logl=25,fig=f3)#,scaletop=0.8)
        plt.draw()
        plt.savefig(prefix+'P_ELD.png',format='png')
        plt.close()
        
        f3=plt.figure()#figsize=(8,12))
        self.plot_2Dfield(component='TDPwEIBW',logl=25,fig=f3)#,scaletop=0.8)
        plt.draw()
        plt.savefig(prefix+'P_IBW.png',format='png')
        plt.close()
        
        for ispec in range(self.nspec):
            f3=plt.figure()#figsize=(8,12))
            self.plot_2Dfield(component='TDPwIF',species=ispec+1,logl=25,fig=f3)
            plt.draw()
            plt.savefig(prefix+'P_IF'+str(ispec+1)+'.png',format='png')
            plt.close()
            
            f3=plt.figure()#figsize=(8,12))
            self.plot_2Dfield(component='TdPwIH',species=ispec+1,logl=25,fig=f3)
            plt.draw()
            plt.savefig(prefix+'P_IH'+str(ispec+1)+'.png',format='png')
            plt.close()
        return


    def read_equigs(self, equigsfile='equigs.data'):
        "Read the equilibrium file created by toric in toricmode='equil',isol=0."
        if self.idebug:
            print ("Using ", equigsfile)

        self.equigs = read_equigsfile(equigsfile)
        return
    


####main block
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import toric_tools
    import sys
    import getopt


# get file name if provided
    iprefix=""
    ifile="TORICLH.cdf"
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hp:f:",["help","prefix=",
                                                          "file="])
    except getopt.GetoptError:
        print ("Accepted flags are help and prefix=")
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-p","--prefix"):
            iprefix=arg
        elif opt in ("-f","--file"):
            ifile=arg
        elif opt in ("-h","--help"):
            print(toric_tools.toric_analysis.__doc__)
            print('run \"help toric_tools.toric_analysis\" for help on whole class')

#Load a run
    LHRun=toric_tools.toric_analysis(toric_name=ifile)
    LHRun.threeplots(prefix=iprefix)


#make sequence of plots, ala the old toric idl driver as an option.


#From ips-wrappers
class toric_file():
    def __init__(self,toric_name='fort.9',lean=True):
        from scipy.io import netcdf_file
        #open cql3d netcdf
        #------------------------------------------------------------------------
        try:
            toric_nc = netcdf_file(toric_name,'r')
        except:
            print('toric_file initialization failed: could not find ncdf: ',
                  toric_name)
            raise Exception(
                'toric_file initialization failed: could not find ncdf: ',
                toric_name)

        #read in cdf dimensions
        #------------------------------------------------------------------------
        self.n_of_field_comp  = np.copy(toric_nc.dimensions['n_of_field_comp'])
        self.n_of_pol_modes   = np.copy(toric_nc.dimensions['n_of_pol_modes'])
        self.n_of_pol_pts     = np.copy(toric_nc.dimensions['n_of_pol_pts'])
        self.n_of_pol_modvac  = np.copy(toric_nc.dimensions['n_of_pol_modvac'])
        self.n_of_rad_elem    = np.copy(toric_nc.dimensions['n_of_rad_elem'])
        self.n_of_ionspec     = np.copy(toric_nc.dimensions['n_of_ionspec'])
        self.dim2_of_celem    = np.copy(toric_nc.dimensions['dim2_of_celem'])
        self.n_of_rad_pts     = np.copy(toric_nc.dimensions['n_of_rad_pts'])
        self.n_of_mhd_modes   = np.copy(toric_nc.dimensions['n_of_mhd_modes'])
        self.n_of_mhd_rad_pts = np.copy(toric_nc.dimensions['n_of_mhd_rad_pts'])
        self.data_plawall     = np.copy(toric_nc.dimensions['data_plawall'])

        #transcripe some dims to shorthand
        self.ndims = self.n_of_field_comp
        self.nmod  = self.n_of_pol_modes
        self.ntt   = self.n_of_pol_pts
        self.nelm  = self.n_of_rad_pts
        self.nmhd  = self.n_of_mhd_rad_pts

        #read in cdf variables
        #------------------------------------------------------------------------
        if(lean==False):
            self.enhcol = np.copy(toric_nc.variables['enhcol'].data)
            self.dnures = np.copy(toric_nc.variables['dnures'].data)
            self.tnures = np.copy(toric_nc.variables['tnures'].data)
            self.ant_length = np.copy(toric_nc.variables['ant_length'].data)
            self.ant_constant = np.copy(toric_nc.variables['ant_constant'].data)
            self.ant_position = np.copy(toric_nc.variables['ant_position'].data)
            self.torus_radius = np.copy(toric_nc.variables['torus_radius'].data)
            self.axis_radius = np.copy(toric_nc.variables['axis_radius'].data)
            self.plasma_radius = np.copy(toric_nc.variables['plasma_radius'].data)
            self.sep_radius = np.copy(toric_nc.variables['sep_radius'].data)
            self.farshield_radius = np.copy(toric_nc.variables['farshield_radius'].data)
            self.antenna_radius = np.copy(toric_nc.variables['antenna_radius'].data)
            self.wall_radius = np.copy(toric_nc.variables['wall_radius'].data)
            self.b_zero = np.copy(toric_nc.variables['b_zero'].data)
            self.b_axis = np.copy(toric_nc.variables['b_axis'].data)
            self.tor_current = np.copy(toric_nc.variables['tor_current'].data)
            self.psi_edge = np.copy(toric_nc.variables['psi_edge'].data)
            self.ppjte = np.copy(toric_nc.variables['ppjte'].data)
            self.ppjti = np.copy(toric_nc.variables['ppjti'].data)
            self.Shfr_shift_axis = np.copy(toric_nc.variables['Shfr_shift_axis'].data)
            self.Shfr_shift_wall = np.copy(toric_nc.variables['Shfr_shift_wall'].data)
            self.ellip_axis = np.copy(toric_nc.variables['ellip_axis'].data)
            self.ellip_edge = np.copy(toric_nc.variables['ellip_edge'].data)
            self.ellip_wall = np.copy(toric_nc.variables['ellip_wall'].data)
            self.triang_edge = np.copy(toric_nc.variables['triang_edge'].data)
            self.triang_wall = np.copy(toric_nc.variables['triang_wall'].data)
            self.vert_shift_axis = np.copy(toric_nc.variables['vert_shift_axis'].data)
            self.vert_shift_wall = np.copy(toric_nc.variables['vert_shift_wall'].data)
            self.vert_triang_edge = np.copy(toric_nc.variables['vert_triang_edge'].data)
            self.vert_triang_wall = np.copy(toric_nc.variables['vert_triang_wall'].data)
            self.edge_skewdness = np.copy(toric_nc.variables['edge_skewdness'].data)
            self.dist_plafars = np.copy(toric_nc.variables['dist_plafars'].data)
            self.dist_plaant = np.copy(toric_nc.variables['dist_plaant'].data)
            self.dist_plawall = np.copy(toric_nc.variables['dist_plawall'].data)
            self.centr_elec_dens = np.copy(toric_nc.variables['centr_elec_dens'].data)
            self.centr_elec_temp = np.copy(toric_nc.variables['centr_elec_temp'].data)
            self.sep_elec_dens = np.copy(toric_nc.variables['sep_elec_dens'].data)
            self.sep_elec_temp = np.copy(toric_nc.variables['sep_elec_temp'].data)
            self.centr_ion_temp = np.copy(toric_nc.variables['centr_ion_temp'].data)
            self.sep_ion_temp = np.copy(toric_nc.variables['sep_ion_temp'].data)
            self.so_thickness = np.copy(toric_nc.variables['so_thickness'].data)
            self.so_dens_length = np.copy(toric_nc.variables['so_dens_length'].data)
            self.so_ele_temp_length = np.copy(toric_nc.variables['so_ele_temp_length'].data)
            self.so_ion_temp_len = np.copy(toric_nc.variables['so_ion_temp_len'].data)
            self.ppnei = np.copy(toric_nc.variables['ppnei'].data)
            self.ppnee = np.copy(toric_nc.variables['ppnee'].data)
            self.pptei = np.copy(toric_nc.variables['pptei'].data)
            self.pptee = np.copy(toric_nc.variables['pptee'].data)
            self.pptii = np.copy(toric_nc.variables['pptii'].data)
            self.pptie = np.copy(toric_nc.variables['pptie'].data)
            self.relef = np.copy(toric_nc.variables['relef'].data)
            self.ielef = np.copy(toric_nc.variables['ielef'].data)
            self.poynt_flux = np.copy(toric_nc.variables['poynt_flux'].data)
            self.pow_prof_elec_fw = np.copy(toric_nc.variables['pow_prof_elec_fw'].data)
            self.pow_prof_eld = np.copy(toric_nc.variables['pow_prof_eld'].data)
            self.pow_prof_ttmpe = np.copy(toric_nc.variables['pow_prof_ttmpe'].data)
            self.pow_prof_mxde = np.copy(toric_nc.variables['pow_prof_mxde'].data)
            self.pow_prof_ibwe = np.copy(toric_nc.variables['pow_prof_ibwe'].data)
            self.pow_prof_tot_elec = np.copy(toric_nc.variables['pow_prof_tot_elec'].data)
            self.pow_prof_ICfund_ions = np.copy(toric_nc.variables['pow_prof_ICfund_ions'].data)
            self.pow_prof_ICharm_ions = np.copy(toric_nc.variables['pow_prof_ICharm_ions'].data)
            self.specif_volume = np.copy(toric_nc.variables['specif_volume'].data)
            self.specif_area = np.copy(toric_nc.variables['specif_area'].data)
            self.tot_rf_current = np.copy(toric_nc.variables['tot_rf_current'].data)
            self.zeff = np.copy(toric_nc.variables['zeff'].data)
            self.prof_rf_curr = np.copy(toric_nc.variables['prof_rf_curr'].data)
            self.equil_file = np.copy(toric_nc.variables['equil_file'].data)
            self.profnt_file = np.copy(toric_nc.variables['profnt_file'].data)

        #small vars needed for stuff like power rescaling, etc.
        self.frequency = np.copy(toric_nc.variables['frequency'].data)
        self.atomic_masses = np.copy(toric_nc.variables['atomic_masses'].data)
        self.atomic_charges = np.copy(toric_nc.variables['atomic_charges'].data)
        self.concentrations = np.copy(toric_nc.variables['concentrations'].data)
        self.psi_mesh = np.copy(toric_nc.variables['psi_mesh'].data)
        self.total_power = np.copy(toric_nc.variables['total_power'].data)
        self.tot_power_eles = np.copy(toric_nc.variables['tot_power_eles'].data)
        self.fw_power_elec = np.copy(toric_nc.variables['fw_power_elec'].data)
        self.ib_power_elec = np.copy(toric_nc.variables['ib_power_elec'].data)
        self.ICfund_power_ions = np.copy(toric_nc.variables['ICfund_power_ions'].data)
        self.ICharm_power_ions = np.copy(toric_nc.variables['ICharm_power_ions'].data)
