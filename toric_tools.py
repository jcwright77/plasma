#!/usr/bin/env python
##!/opt/bin/python
#upgrading to scipy version 0.100dev , changing interface
import numpy as np
import numpy.fft as ft
import scipy.io.netcdf as nc
import matplotlib.pyplot as plt
import os
from matplotlib import ticker, cm
#other deps below
    #import f90nml
    #from periodictable import elements    

# Plasma species in template namelist
def get_spec_toric(toricnml):
    "Collect info on species for ICRF sim in TORIC in nice readable format"
    #import f90nml
    from periodictable import elements    
    
    spec_toric=list(zip(map(round,toricnml['equidata']['atm']),map(int,toricnml['equidata']['azi'])))
    for i,s in enumerate(spec_toric):
      name=str(elements[s[1]][s[0]])
      spec_toric[i]={'name':name,'A':spec_toric[i][0],'Z':spec_toric[i][1],'Conc%':100*toricnml['equidata']['aconc'][i]}
    
    spec_toric.insert(0,{'name':'e', 'A':0, 'Z':-1 , 'Conc%': 100})
    return spec_toric


def print_vector(nrep,fstr,a):
    """
    Converts an array of numbers into a string formated by fstr with nrep values per line.
    """
    n=a.size
    pa=""
    ta=a.reshape((a.size),order='F')
    for j in range(0,n,nrep):
        pa=pa+ "".join(map(lambda f: fstr % f, ta[j:min(j+nrep,n)]))+"\n"
    return pa


def write_profnt(namelist,equidt,version='profnt2'):
    """
    Inputs:
        namelist: as created by f90nml from toric.inp
        equidt: python dictionary of profiles to write out
            psipro: sqrt(Psipol/Psipol[a])
            tbne: [cm-3] on psipro mesh
            tbte: [keV]  on psipro mesh
            iatm: array of atomic masses/ (C12/12)
            iazi: array of atomic numbers
            tbni: ion densities on psipro mesh
            tbi_provv: ion temperatures on psipro mesh
            nspec: number of ion species
            
         version: Generally format2 is used. File is self describing in
                  number of elements and profiles.
    """

    filename=namelist['equidata']['profnt_file']
    with open(filename,'w') as of:
        nprodt=equidt['psipro'].size
        profiles=['psipro','tbne','tbte','tbi_provv']
        if version=='profnt1':
            for profile in profiles:
                of.write('{:<10s}{:4d}\n'.format(profile, nprodt ))
                of.write(print_vector(5,'%16.9e',equidt[profile]))

        if version=='profnt2':
            nspec=len(equidt['iatm'])
            mainsp=1
            namelist['equidata']['mainsp']=mainsp
            kdiff_idens=equidt['kdiff_idens'] #0 #specify concentrations
            kdiff_itemp=equidt['kdiff_itemp'] #0 #one ion temp
            of.write('{:<10s}{:4d}{:4d}{:4d}{:4d}{:4d}\n'.format('profnt_py',nprodt,nspec,
                                                                 mainsp,kdiff_idens,kdiff_itemp))
            for isp in range(nspec):
                of.write('{:4d}{:4d}\n'.format(int(equidt['iatm'][isp]),int(equidt['iazi'][isp])) )
            profiles=['psipro','tbne','tbte']
            for profile in profiles:
                of.write('{:<10s}{:4d}\n'.format(profile, nprodt ))
                of.write(print_vector(5,'%16.9e',np.array(equidt[profile])))

            #write ion densities and temperatures
            for isp in range(nspec):
                if kdiff_idens==0:
                    of.write('{:<10s}\n'.format('ni_conc'+str(isp)))
                    of.write('%16.9e \n' % equidt['tbni'][isp]) 
                else:
                    of.write('{:<10s}\n'.format('tbni'+str(isp)))
                    of.write(print_vector(5,'%16.9e',equidt['tbni'][:,isp]))

                if kdiff_itemp==0 and isp==0:
                    of.write('{:<10s}\n'.format('ion_temp') )
                    of.write(print_vector(5,'%16.9e',equidt['tbi_provv'])) 

                if kdiff_itemp==1:
                    of.write('{:<10s}\n'.format('ion_temp'+str(isp)) )
                    of.write(print_vector(5,'%16.9e',equidt['tbi_provv'][:,isp])) 

#this routine is still incomplete
def read_equidt(filename):
    import fortranformat as ff
    equidt={}
    with open(filename,'r') as of:
        line = of.readline()
        reader = ff.FortranRecordReader('(A10,5i4)')
        var_name, nprodt, nspec, mainsp,kdiff_idens, kdiff_itemp = reader.read(line)
        
            
        line = of.readline()
        reader = ff.FortranRecordReader('(A10,5i4)')

        if var_name=='Rfxqlo_Pro':
        #  DMC -- define flag to indicate TRANSP format variant
            print(' Detected: TRANSP "Rfxqlo_Pro" file format variant!')
            kdiff_itemp = 1
            kdiff_idens = 1
            rfxqlo_pro_variant = True
        else:
            rfxqlo_pro_variant = False

        if kdiff_itemp==0:
            nsptmp = 1
        else:
            nsptmp = 10 # place holder nspec



def formattedwrite(file,a):
  sza=len(a)
  for idx in range(0,int(sza/4)*4,4):
    file.write(f"{a[idx]:18.9E}{a[idx+1]:18.9E}{a[idx+2]:18.9E}{a[idx+3]:18.9E}\n")
  rem = np.mod(sza,4)
  if rem>0:
    file.write(''.join([ "%18.9E" % x for x in a[-rem:] ])+"\n") #print remaining elements

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
  with open(equigsfile, "w") as file:
      torlheq_mmodes=20
      torlheq_psimodes=eq['xmap'].shape[0]
      psimap=eq['psipolmap']
      equigs={}

      file.write(' Major radius (central)(m) [generated by plasma.equigs.py]\n')
      file.write(f"{eq['rcentr']:18.9E}\n")
      equigs["rtorm"] =eq['rcentr']

      file.write(' Major radius (axis)(m)\n')
      file.write(f"{eq['rmaxis']:18.9E}\n")
      equigs["rmaxis"] =eq['rmaxis']

      file.write(' Magnetic field at major radius (m)\n')
      file.write(f"{eq['bcentr']:18.9E}\n")
      equigs["bcentr"] =eq['bcentr']

      file.write(' Total toroidal current (kA)\n')
      equigs["torcur"] =eq['current']/1000.
      file.write(f"{equigs['torcur']:18.9E}\n") #eqdsk is in Amps, torlh in kAmps
      
      file.write(' Number of poloidal modes\n')
      equigs["imom"] = torlheq_mmodes
      file.write(f"{torlheq_mmodes:5}\n")

      file.write(' Number of radial mesh points\n')
      equigs["nmhd"] = torlheq_psimodes
      file.write(f"{torlheq_psimodes:5}\n")

      file.write(' Radial mesh\n')
      rhopol=eq['rhopolmap']
      equigs["srad"] = rhopol
      formattedwrite(file,rhopol)

      file.write(' Fourier equilibrium coefficients\n')
      equigs["rzmcs2d"]=eq['rzmc2d']
      rmc2d,rms2d,zmc2d,zms2d=eq['rzmc2d']
      formattedwrite(file,rmc2d[:,0]) #dc modes
      formattedwrite(file,zmc2d[:,0])
      for i in range(1,rmc2d.shape[1]):
        formattedwrite(file,rmc2d[:,i])
        formattedwrite(file,zms2d[:,i])
        formattedwrite(file,rms2d[:,i])
        formattedwrite(file,zmc2d[:,i])
  
      file.write(' Safety factor\n')
      qmap=np.interp(eq['psipolmap'],eq['fluxGrid'],eq['qpsi'])
      equigs["qqf"] = qmap
      formattedwrite(file,qmap)

      file.write(' Current profile [kA]\n')
      equigs["Icurr"]=eq['Ipsimap']/1000.0
      formattedwrite(file,equigs["Icurr"])

      file.write(' Covariant B_phi, R*B_phi (m*T)\n')
      gmap=np.interp(eq['psipolmap'],eq['fluxGrid'],eq['fpol'])
      equigs["gcov"] = gmap
      formattedwrite(file,gmap)

      file.write(' Rho toroidal\n')
      equigs["rhotor"]=eq['rhotormap']
      formattedwrite(file,eq['rhotormap'])

      file.write(' Fraction Psi poloidal at last surface\n')
      equigs["lastpsi"]=eq['lastpsi']
      file.write(f"{equigs['lastpsi']:18.9E}\n")

  return equigs




class toric_analysis:
    """Class to encapsulate tools used for toric analysis.
    Works with scipy 0.8.0 and scipy 0.10.0 and numpy 1.6

    Typical invocation:
    import toric_tools
    R=toric_tools.toric_analysis('toric.ncdf',mode='ICRF') #'LH' for lower hybrid
    R.plot_2Dfield(component='Re2Ezeta',logl=10) #etc

    Important functions:
    R.info() # netcdf contents and metadata
    R.plotpower(power='PwIF',species=1) #different power profiles
    R.plot_2Dfield(component='Re2Ezeta',logl=10) #Two Dim plots of quantities
    R.plot_1Dfield(component='Re2Ezeta') #One Dim plots of quantities
    R.threeplots() #2D Ez, electron power and poynting flux and poloidal spectrum, works for LH only, presently
    """


    def __init__ (self, toric_name='None', mode='LH',
        idebug=False, comment='', layout='poster', path="./"):
        import socket
        from time import gmtime

        __version__=1.0

        self.toric_name=toric_name
        self.mode = mode
        self.__version__ = __version__
        self.idebug = idebug

        self.mylw=1.0
        self.mypt=18.0
        self.fsc=2.0
        self.fw='bold'
        self.set_layout(layout)

        self.path=path

        self.prov = {"user":"noname","host":"noname","gmtime":"notime","runid":"noid",
                "path":"", "comment":""}
        self.label = True
        self.equigs = {}
        self.toricdict={}

        if (self.mode[:2]=='LH'):
            self.namemap={'xpsi':'tpsi','poynt':'vpoynt','pelec':'S_eld',
                         'e2d_z':'E2d_z_re','xplasma':'x_plasma', 'zplasma':'z_plasma',
                          'xeqpl':'xeqpl'}
            if self.toric_name=='None': self.toric_name='TORICLH.cdf'
        else:
            self.namemap={'xpsi':'Pw_abscissa','poynt':'PoyFlx','pelec':'PwE',
                         'e2d_z':'Re2Ezeta','xplasma':'Xplasma', 'zplasma':'Zplasma',
                          'xeqpl':'Ef_abscissa'}
            if self.toric_name=='None': self.toric_name='toric.ncdf'
##Open the toric netcdf file read only
        try:
            self.cdf_hdl = nc.netcdf_file(path+self.toric_name,mmap=False )#,'r')
        except IOError:
            print ('CRITICAL: ',self.toric_name,' not found.')
            self.cdf_hdl = -1

        try:
            self.qlde_hdl = nc.netcdf_file(path+"toric_qlde.cdf",mmap=False)#,'r')
        except IOError:
            print ('Non-CRITICAL: ',path+"toric_qlde.cdf",' not found.')
            self.qlde_hdl = -1

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
            
        return

    
    def info( self ):
        "Prints a list of the contents of present TORIC3D output files"

        if (self.cdf_hdl != -1):
            for hdl in [self.cdf_hdl]:
                print ('The toric file, ',self.toric_name,', contains:')
                print ('----------------------------------------------')
                print ("The global attributes: ",hdl.dimensions.keys())        
                print ("File contains the variables: ", hdl.variables.keys())

        if (self.qlde_hdl != -1):
            for hdl in [self.qlde_hdl]:
                print ('The toric file, ',self.toric_name,', contains:')
                print ('----------------------------------------------')
                print ("The global attributes: ",hdl.dimensions.keys()  ) 
                print ("File contains the variables: ", hdl.variables.keys())


        print ('----------------------------------------------')
        print ("Provenance metadata: ", self.prov)

        return


    def toricparam( self ):
        "Fill a dictionary of toric scalar values"

        self.toricdict['comment']='Dictionary of relevant toric scalars'
        self.toricdict['mode']=self.mode
        self.toricdict['nelm']=self.cdf_hdl.dimensions['nelm']
        self.toricdict['nptpsi']=self.cdf_hdl.dimensions['nptpsi']
        self.toricdict['ntt']=self.cdf_hdl.dimensions['ntt']

        return


    def plotb0( self, ir=45, db=0,eps=0 ):
        """Contour plot of bounce averaged Dql coefficient, dB0 """
        if (self.qlde_hdl == -1):
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
        return
            

    def plotpower( self, xaxis=None, power=None, species=None ):
        """Plot power profiles versus specified radius for all, or listed species.
           Overplots by default. species is 0 indexed.
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
            plt.xlabel(r'$\sqrt{\psi_{pol}}$')
        cf=plt.gcf()
        cf.subplots_adjust(bottom=0.14)

        return l

    
    def psiplot( self, y ):
        "Plot versus rhopsi. Returns handle on line to modify line style if desired using setp."
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
        """Mapping toric field in celef.cdf to re and im parts of the three components
        for LH mode.
        """

        return


    def __getvar__( self, name ):
        """Internal function to retrieve variable from data file with checking.
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

    def spectrum( self, component='undef',maxr=1.,cx=0,levels=-1, q=None ):
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
        if (np.size(levels)==1):
            nlevels=7
#            levels=np.arange(nelm/nlevels,nelm-1,nelm/nlevels)
            levels=(np.arange(nlevels)*nelm/nlevels).astype(int)
        else:
            levels=(np.array(levels)*nelm).astype(int)
            nlevels=np.size(levels)

        levels=levels[1:nlevels]
#        levels=nelm-np.arange(1,20,2)
        rlevels=rad[levels]

        th = np.arange(ntt)-ntt/2

        ymax = 0.
        ymin = 0.

        i=0
        thq=th
        for indr in range(levels.size): #levels:  #fft in python isn't normalized to N
            ir=levels[indr]
            #print ('levels',ir,levels[indr],rlevels[indr],th)
            if q!=None:
                thq=-2.5*(1+0.3)/(1+0.3*rlevels[indr])*(1+th/191./q(rlevels[indr]))

            #print (thq)
            ffield = ft.fftshift(np.log10(abs(ft.fft(field[:,ir]))/float(ntt)+1.e-20))
            ymax = np.max( [ymax, np.max(ffield)] )
            ymin = np.min( [ymin, np.min(ffield)] )
            plabel='%5.2f' % rad[ir]
            if self.bw:
                plt.plot( thq, ffield, label=plabel, linestyle=self.ls[i],color='k')
                i=i+1
            else:
                plt.plot( thq, ffield, label=plabel )

        ffield = ft.fftshift(np.log10(abs(ft.fft(field[:,nelm-1]))/float(ntt)+1.e-20))
        ymax = np.max( [ymax, np.max(ffield)] )
        ymin = np.min( [ymin, np.min(ffield)] )
#plot antenna spectrum
        plabel='ant'
        print ("range, levels", rlevels)
        print ("ymax", ymax,ymin)
        plt.plot( thq, ffield,  label=plabel, color='grey',linewidth=2 )
        cf=plt.gcf()
        cf.subplots_adjust(right=0.76)
        plt.axis ('tight')
        if q!=None:
            plt.axis( xmin=-8,xmax=8 )
        else:
            plt.axis( xmin=-ntt/4, xmax=ntt/4)
        plt.axis( ymin=-10)
        plt.legend(loc=(1.05,0))
        plt.xlabel('m')
        plt.ylabel('log10 scale')
        plt.title('Poloidal spectrum on labeled flux surfaces')
        plt.draw()
        return

    def set_layout( self, layout='poster' ):

        if (layout == 'paper'):
            self.mylw=2.0
            self.mypt=10.0
            self.fsc=1.0
            self.fw='normal'

        params = {
            'axes.linewidth': self.mylw,
            'lines.linewidth': self.mylw,
            'axes.labelsize': self.mypt,
            'font.size': self.mypt,
            'legend.fontsize': self.mypt,
            'axes.titlesize': self.mypt+2.0,
            'xtick.labelsize':self.mypt,
            'ytick.labelsize':self.mypt,
            'font.weight'  : self.fw,
            'text.usetex' : False
            }
        plt.rcParams.update(params)

        return


#note that if plot commands are in the toplevel, they will not return
#to the prompt, but wait to be killed.
    def plot_2Dfield(self, component='E2d_z',species=None,logl=0,xunits=1.0,axis=(0.0,0.0),
                     im=False, scaletop=1.0, scalebot=1.0,ax='undef',fig='undef',
                     maxsurface=0.99,lscaletop=0.0,lscalebot=0.0):
        """
    
        example of using netcdf python modules to plot toric solutions
        requires numpy and matplotlib and netcdf modules for python.

        Note, under windows you need netcdf.dll installed in SYSTEM32 and the file
        system cannot follow symbolic links.  The DLL needs to have executable
        permissions.

        To overplot with limiter, made from efit plotter:
        R.plot_2Dfield(component='Im2Eplus',logl=20,xunits=0.01,axis=maxis,fig=fig1)

        Easier is to plot solution first, then overplot limiter, scaled appropriately:
        p.plot ( rlim*100.-maxis[0], zlim*100.-maxis[1], 'k', linewidth = 2 )
        
        """

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

        print (self.mode,'mode')
        if (self.mode[:2]=='LH'):
            if (im):
                im_e2dname=component+'_im'
                title='|'+component+'|'
                component=component+'_re'
        else:
            if (component=='E2d_z'):
                component='Ezeta'

            im_e2dname='Im2'+component
            if (im) : title='|'+component+'|'
            component='Re2'+component


#note change to use ().data instead of np.array() in scipy0.8.0
        if (component=="power" and self.mode[:2]=='LH'):
            e2d = self.get_power2D()
        else:
            e2d = (self.cdf_hdl.variables[component]).data

        if (im):
            im_e2d=(self.cdf_hdl.variables[im_e2dname]).data 
            e2d = abs(e2d+1.0j*im_e2d)

        if (self.mode[:2]!='LH' and species):
            print('plot2D, indexing species', species)
            e2d = e2d[:,:,species-1]


        print ("2D Matrix shape:", np.shape(xx))


    #contour with 3 args is confused unless arrays are indexed slices
    #need wrapper to close periodicity in theta direction for this
    #tricky, array indexing different from ncdf slicing
    #this step is needed because periodic dimension is not closed.
    #i.e. its [0,pi) not [0,pi]
        dd=np.shape(xx)
        sx=dd[0]
        sy=dd[1]
        lastpsi=int(sy*maxsurface); print(sx,sy,lastpsi)
        
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
        
        emax=np.max(ee2d[:,:lastpsi].ravel()) #[ee2d.argmax()]
        emin=np.min(ee2d[:,:lastpsi].ravel())#[ee2d.argmin()]

    #contouring levels
        rmax=max([abs(emax),abs(emin)])*scaletop
        rmin=min([0.,emax,emin])*scalebot
        #val=arange(emin,emax,(emax-emin)/25.,'d')
        val=np.arange(-rmax*1.1,rmax*1.1,(rmax+rmax)/25.,'d')
        if (im):
            val=np.arange(rmin,rmax*1.1,(rmax)/24.,'d')
            print ("values",val)

    #reverse redblue map so red is positive
           # revRBmap=cmap_xmap(lambda x: 1.-x, cm.get_cmap('RdBu'))

    #finally, make the plot
        cwidth=xxx.max()-xxx.min()
        cheight=yyy.max()-yyy.min()
        asp=cheight/cwidth
        print ("plot aspect ratio:", asp)

 #leave space for bar
        if (fig=='undef'):
            fig=plt.figure(figsize=(self.fsc*3.0+legend_frac,3.0*self.fsc*asp))
            fig.subplots_adjust(left=0.02,bottom=0.15,top=0.90)

        sax=plt.axes().set_aspect(1, 'box') # the right way to control aspect ratio

        maxpsi=xxx.shape[1]-1
        plt.plot(xxx[:,maxpsi],yyy[:,maxpsi],'k-')

        #add LCF
        lcfpsi=self.cdf_hdl.dimensions['PsiPwdDim']
        plt.plot(xxx[:,lcfpsi],yyy[:,lcfpsi],'grey')
        
    #read ant length.  Calculate arc length vs theta to this value/2
    #in each direction, this plots the antenna location
        anthw=max(int(sx*0.01),4)
        plt.plot(xxx[sx-anthw+1:sx+1,maxpsi],yyy[sx-anthw+1:sx+1,maxpsi],'g-',linewidth=6)
        plt.plot(xxx[0:anthw,maxpsi],yyy[0:anthw,maxpsi],'g-',linewidth=6)
        print("antenna: ", yyy[sx-anthw+1:sx+1,maxpsi])
        if self.label:
            ax=plt.gca()
            sublabel=self.prov['path']
            print (sublabel)
            plt.text(-0.2,-0.2,sublabel,transform = ax.transAxes)

        print ("interactive off while plotting")
#        plt.ioff()

        if (logl > 0):
            title='log10 '+title
            barfmt='%3.1f'

    ##labels and titles
    #xlabel(getattr(xx,'long_name')+'('+getattr(xx,'units')+')')
    #ylabel(getattr(yy,'long_name')+'('+getattr(yy,'units')+')')
    #title(getattr(e2d,'long_name')+'('+getattr(e2d,'units')+')')
        plt.xlabel('X(cm)')
        plt.ylabel('Z(cm)')
#        plt.title(r'$Re E_{||}$',fontsize=self.mypt+2.0)
        plt.title(title,fontsize=self.mypt+2.0)


        if (logl <= 0):
            CS=plt.contourf(xxx,yyy,ee2d,val,cmap=cm.jet) #30APR2009 removed *0.2

        if (logl > 0):
#            lee2d=np.sign(ee2d)*np.log(np.sqrt(np.abs(ee2d)**2+1)+np.abs(ee2d))/np.log(10)
            lee2d=np.log(np.abs(ee2d)+1.0)/np.log(10)
            rmax=lee2d.ravel()[lee2d.argmax()]+lscaletop
            rmin=lee2d.ravel()[lee2d.argmin()]+lscalebot
            val=np.arange(rmin,rmax,(rmax-rmin)/(logl*1.0),'d')
            CS=plt.contourf(xxx,yyy,lee2d,val,cmap=cm.jet)
        print ("interactive on")
#        plt.ion()
##put the contour scales on the plot
#tricky, fraction needs to be specified to be part by which horizontal exceed vertical

        cbar=plt.colorbar(CS,format=barfmt,ax=sax)
        cbar.ax.set_ylabel('levels')

        print ("contour values",CS.levels,'xx',rmax,rmin)

        return CS,cbar

  
        
### user routines using the above, could be in a different module
    def powpoynt( self ):
        "Plots powers and poynting flux"
        fig = plt.figure(figsize=(16,9) )
        ax1 = fig.add_subplot(111)
        line1,=self.psiplot(self.namemap['pelec'])
#can use setp(lines, ) to change plot properties.
        plt.setp(line1,color='b' ,label=r'$P_{eld}$')
     #   line1.setlabel('<ExB>')


#add first two species if ICRF, add logic to plot if power percent is larger than 0.5%
        if (self.mode[:2]!='LH'):
           line2,=self.plotpower(power='PwIF',species=1)
#           line3,=self.plotpower(power='PwIF',species=2)
         
           line4,=self.plotpower(power='PwIF',species=3)
           plt.setp(line4,color='g',label='Fund sp3') #add species name
#           line5,=self.plotpower(power='PwIH',species=1)
           line6,=self.plotpower(power='PwIH',species=2)
           plt.setp(line6,color='k',label='Harm sp2') #add species name
#           line7,=self.plotpower(power='PwIH',species=3)

        ax2 = ax1.twinx()
        line2,=self.psiplot(self.namemap['poynt'])
#set axis floor at 0
        #plt.gca().set_ylim(0)
        ax1.set_ylim(0)
        ax2.set_ylim(0)
        ax1.set_ylabel('Power',color='b')        
#change color and symbol
        plt.setp(line2,color='r', label='<ExB>')
        ax2.set_ylabel('Poynting',color='r')
        ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

#make  legend too
        plt.legend( handles=[line1,line2,line4,line6], loc='center right', 
                     ncol=1, fancybox=True, shadow=True)
        
        plt.tight_layout()
        plt.draw()
        return fig

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
        fig.subplots_adjust(left=0.12,bottom=0.12,top=0.96,right=0.82,hspace=0.32)
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
#            cmd="sed -n '/elec/,/,/ s/-0\./ -0./pg'  < torica.sol > torica_2dpower.sol"
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
        """Makes and saves the three most commonly used plots. Plots are saved in 
        the current directory. An optional prefix can be used to label them or change
        the save path.
        * Power and poynting flux on one plot as eps.
        * The polodial power spectrum on six flux surfaces for convergence as eps.
        * And the 2D parallel electric field contour plot as a png."""

        self.spectrum(cx=1)
        plt.draw()
        plt.savefig(prefix+'spectrum.pdf',format='pdf')
        plt.savefig(prefix+'spectrum.png',format='png')

        if (self.mode[:2]!='LH'):        
            self.plot_2Dfield(component='Eplus', maxsurface=0.9,im=True,logl=25)
            plt.draw()
            plt.savefig('log10Eplus2d.png',format='png')
            self.plot_2Dfield(component='Eplus',maxsurface=0.9)#,scaletop=.8)
            plt.draw()
            plt.savefig('Eplus2d.png',format='png')
            
        
        self.plot_2Dfield(im=True,logl=25)#,scaletop=0.8)
        plt.draw()
        plt.savefig(prefix+'log10Ez2d.png',format='png')

        self.powpoynt()
        plt.draw()
        plt.savefig(prefix+'powerpoynt.pdf',format='pdf')
        plt.savefig(prefix+'powerpoynt.png',format='png')

        return


    
#Handling equigs file
    def __get_varname(self, f):
        "Reads next line from file f and returns it, optionally printing it."
        varname=f.readline()
        if self.idebug:
            print (f.name,varname)
        return varname

    
    def read_equigs(self, equigsfile='equigs.data'):
        "Read the equilibrium file created by toric in toricmode='equil',isol=0."
        if self.idebug:
            print ("Using ", equigsfile)
        equigs_hdl=open(equigsfile,'r')

        varname = self.__get_varname(equigs_hdl)
        self.equigs["rtorm"] = np.fromfile(equigs_hdl,sep=" ",count=1,dtype=float)[0]

        varname = self.__get_varname(equigs_hdl)
        self.equigs["raxis"]= np.fromfile(equigs_hdl,sep=" ",count=1,dtype=float)[0]

        varname = self.__get_varname(equigs_hdl)
        self.equigs["bzero"] = np.fromfile(equigs_hdl,sep=" ",count=1,dtype=float)[0]

        varname = self.__get_varname(equigs_hdl)
        self.equigs["torcur"]= np.fromfile(equigs_hdl,sep=" ",count=1,dtype=float)[0]

        varname = self.__get_varname(equigs_hdl)
        self.equigs["imom"] = np.fromfile(equigs_hdl,sep=" ",count=1,dtype=int)[0]
        imom = self.equigs["imom"]

        varname = self.__get_varname(equigs_hdl)
        self.equigs["nmhd"] = np.fromfile(equigs_hdl,sep=" ",count=1,dtype=int)[0]
        nmhd = self.equigs["nmhd"]

        varname = self.__get_varname(equigs_hdl)
        self.equigs["srad"] = np.fromfile(equigs_hdl,sep=" ",count=nmhd,dtype=float)

        #this needs to be reshaped or remapped into the R,Z sin cos arrays toric uses
        varname = self.__get_varname(equigs_hdl)
        self.equigs["rzmcs2d"] = np.fromfile(equigs_hdl,sep=" ",
                                                count=2*nmhd+4*nmhd*imom,dtype=float)

        varname = self.__get_varname(equigs_hdl)
        self.equigs["qqf"] = np.fromfile(equigs_hdl,sep=" ",count=nmhd,dtype=float)

        #logic checking for "END"
        varname = self.__get_varname(equigs_hdl)
        self.equigs["jcurr"] = np.fromfile(equigs_hdl,sep=" ",count=nmhd,dtype=float)

        varname = self.__get_varname(equigs_hdl)
        self.equigs["gcov"]= np.fromfile(equigs_hdl,sep=" ",count=nmhd,dtype=float)

        varname = self.__get_varname(equigs_hdl)
        self.equigs["rhotor"] = np.fromfile(equigs_hdl,sep=" ",count=nmhd,dtype=float)

        varname = self.__get_varname(equigs_hdl)
        self.equigs["lastpsi"] = np.fromfile(equigs_hdl,sep=" ",count=1,dtype=float)[0]

        equigs_hdl.close()


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
        opts, args = getopt.getopt(sys.argv[1:], "hp:f:",["help","prefix=","file="])
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



