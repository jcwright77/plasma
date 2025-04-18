#!/usr/bin/env python
#
#   Try reading g-eqdsk file with re (regex) module
#   instead of the non-existant fortran format code
#   python feature.
#
#   WARNING: this code has only been testing on the 
#   two files listed below and my regex skills are
#   quite poor (as is my python ) so you have been 
#   warned.
#
#   DLG - 14-Aug-08
#
#   JCW - 27-Jan-10
#
# some improvements to regular expression and implemented direct casting
# of string arrays to floating arrays.
# Regexp are a bit fragile here. Should take advantage of the known data
# ordering with n.fromfile
#
#   JCW - 02-May-12 
#
#   Make into a function call returning a structure with optional outputs
#
# some eqdsk files have data in 16.8, some in 16.9

#Courtesy of OMFIT eq_JET.py
# Need following data:
# GEQDSK NAME                                                          - EQUIVALENT EFIT PPF NAME:
# RDIM: Horizontal dimension in meter of computational box             - FROM PSIR - R Grid
# ZDIM: Vertical dimension in meter of computational box               - FROM PSIZ - Z Grid
# RLEFT: Minimum R in meter of rectangular computational box           - AS ABOVE
# ZMID: Z of center of computational box in meter                      - AS ABOVE
# RMAXIS: R of magnetic axis in meter                                  - RMAG
# ZMAXIS: Z of magnetic axis in meter                                  - ZMAG
# SIMAG: poloidal flux at magnetic axis in Weber /rad                  - FAXS (W/rad)
# SIBRY: poloidal flux at the plasma boundary in Weber /rad            - FBND (W/rad)
# RCENTR: R in meter of vacuum toroidal magnetic field                 - 2.96m
# BCENTR: Vacuum toroidal magnetic field in Tesla at RCENTR            - BVAC
# CURRENT: Plasma current in Ampere                                    - (XIP = measured) XIPC - Calculated plasma current
# FPOL: Poloidal current function in m-T,  F = RBT  on  flux grid      - F
# PRES: Plasma pressure in nt / m2 on uniform flux grid                - P
# FFPRIM: FF'(PSI) in (mT)2 / (Weber /rad) on uniform flux grid        - DFDP * Mu0
# PPRIME: P'(PSI) in (nt /m2) / (Weber /rad) on uniform flux grid      - DPDP
# PSIZR: Poloidal flux in Weber / rad on the rectangular grid points   - PSI
# QPSI: q values on uniform flux grid from axis to boundary            - Q
# NBBBS: Number of boundary points                                     - NBND
# LIMITR: Number of limiter points                                     - length(RLIM)
# RBBBS: R of boundary points in meter                                 - RBND
# ZBBBS: Z of boundary points in meter                                 - ZBND
# RLIM: R of surrounding limiter contour in meter                      - RLIM
# ZLIM: Z of surrounding limiter contour in meter                      - ZLIM

def readGEQDSK(filename='eqdsk.dat', dointerior=False, doplot=None, width=9, 
           	cocos=1, dolimiter=None, ax=None, dodebug=False):
    import re
    import numpy as np
    import pylab as plt

    file = open (filename)
    data    = file.read ()

    dimensionsRE    = re.compile ( ' {1,3}\\d?\\d?\\d?\\d\\d' ) # Equivilant to i5 fortran code, JCW these should be i4
    dimensionsRE4    = re.compile ( ' {1,3}\\d?\\d?\\d?\\d' ) # Equivilant to i4 fortran code
    headerRE    = re.compile ( '^.*\\n') # First line
    if width==9:
        valuesRE   = re.compile ( '([ \\-]\\d\\.\\d{9}[eEdD][\\+\\-]\\d\\d)' )   # Equivilant to e16.9 fortran code
    else:
        valuesRE   = re.compile ( '([ \\-]\\d\\.\\d{8}[eEdD][\\+\\-]\\d\\d)' )   # Equivilant to e16.8 fortran code

#bbbsRE  = re.compile ( '( {1,3}\\d?\\d?\\d?\\d\\d {1,3}\\d?\\d?\\d?\\d\\d)' )   # Candidate dimension lines (2i5 fortran code)
    bbbsRE  = re.compile ( r'(?m)^.{10}\n' ) #there should be only one 10 character line

    dataStr     = valuesRE.findall ( data )
    headerStr   = headerRE.findall ( data )
    bbbStr  = bbbsRE.findall ( data )

    file.close ()
    if len(bbbStr) > 0:
        nbbbsStr    = dimensionsRE4.findall ( bbbStr[0] )
        nbbbs   = int ( nbbbsStr[-2] )
        limitr   = int( nbbbsStr[-1] )
    else:
        print('no bounding box found or limiter. should be Line with 2 integers length of 10 characters')
        nbbbsStr = []
        nbbbs = 0  #should be there but cont if not
        limitr = 0
        
    nWnHStr = dimensionsRE4.findall ( headerStr[0][48:] )

    nW  = int ( nWnHStr[1] )
    nH  = int ( nWnHStr[2] )
    #idummy used as 1D size: nV, nW, nH
    nV  = int ( nWnHStr[0] )
    if nV <= 0:  nV = nW  #disable non standard and meaningless
    nV=nW
    
    rdim    = float ( dataStr[0] )
    zdim    = float ( dataStr[1] )

    if dodebug: 
        print("Data string header:", dataStr[0:20] ),
        print("nWnStr string header:", nWnHStr,headerStr, len(headerStr))
        print("Dimensions:", nW, nH, nV, nbbbs, limitr, rdim, zdim )
        print("Size of data:", len(dataStr), 20+nV*5+nW*nH+2*nbbbs+2*limitr  )

    rcentr  = float ( dataStr[2] )
    rleft   = float ( dataStr[3] )
    zmid    = float ( dataStr[4] )

    rmaxis  = float ( dataStr[5] )
    zmaxis  = float ( dataStr[6] )
    simag   = float ( dataStr[7] )
    sibry   = float ( dataStr[8] )
    bcentr  = float ( dataStr[9] )

    current = float ( dataStr[10] )

    fpol    = np.zeros ( nV )
    pres    = np.zeros ( nV )
    ffprim  = np.zeros ( nV )
    pprime  = np.zeros ( nV )
    psizr   = np.zeros ( ( nW, nH ) )
    qpsi    = np.zeros ( nV )
    rbbbs   = np.zeros ( nbbbs )
    zbbbs   = np.zeros ( nbbbs )
    rlim    = np.zeros ( limitr )
    zlim    = np.zeros ( limitr )


#   If you know how to cast a list of strings to
#   a numpy array without a loop please let me 
#   know, as these loops should not be required.

#   1D arrays

    for i in np.arange ( nV ) : 
    
        fpol[i] = dataStr[np.asarray(i+20,dtype=int)]
        pres[i] = dataStr[np.asarray(i+20+nV,dtype=int)]
        ffprim[i] = dataStr[np.asarray(i+20+2*nV,dtype=int)]
        pprime[i] = dataStr[np.asarray(i+20+3*nV,dtype=int)]
        qpsi[i] = dataStr[np.asarray(i+20+4*nV+nW*nH,dtype=int)]

    if dodebug: print('one D arrays: ', fpol[-1],pres[-1], ffprim[-1], pprime[-1], qpsi[-1] )
    for i in np.arange ( nbbbs ) :  
        rbbbs[i]    = dataStr[np.asarray(i*2+20+5*nV+nW*nH,dtype=int)]
        zbbbs[i]    = dataStr[np.asarray(i*2+1+20+5*nV+nW*nH,dtype=int)]
  

    for i in np.arange ( limitr ) :
       
        rlim[i] = dataStr[np.asarray(i*2+20+5*nV+nW*nH+2*nbbbs,dtype=int)]
        zlim[i] = dataStr[np.asarray(i*2+1+20+5*nV+nW*nH+2*nbbbs,dtype=int)]

#   2D array

    for i in np.arange ( nW ) :
        for j in np.arange ( nH ) :
            psizr[i,j] = dataStr[np.asarray(i+20+4*nV+j*nW,dtype=int)]

    rStep   = rdim / ( nW - 1 )
    zStep   = zdim / ( nH - 1 )
    fStep   = -( simag - sibry ) / ( nW - 1 )

    r   = np.arange ( nW ) * rStep + rleft
    z   = np.arange ( nH ) * zStep + zmid - zdim / 2.0

    fluxGrid    = np.arange ( nW ) * fStep + simag

#   Find indices of points inside and outside
#   the rbbbs/zbbbs boundary.
    import matplotlib.path as mplPath
    import numpy as np
    lcf=mplPath.Path( np.column_stack( (rbbbs,zbbbs) ) )
    iiInsideA   = np.zeros ( psizr.shape )
    iiInside = -1
    iiOutside = -1
    if (dointerior):
        for i in np.arange ( nW ) :
            for j in np.arange ( nH ) :
                if lcf.contains_point( (r[i],z[i]) ):
                    iiInsideA[i,j] = 1
                #q1  = np.size ( np.where ( ( r[i] - rbbbs > 0 ) & ( z[j] - zbbbs > 0 ) ) )
                #q2  = np.size ( np.where ( ( r[i] - rbbbs > 0 ) & ( z[j] - zbbbs <= 0 ) ) )
                #q3  = np.size ( np.where ( ( r[i] - rbbbs <= 0 ) & ( z[j] - zbbbs > 0 ) ) )
                #q4  = np.size ( np.where ( ( r[i] - rbbbs <= 0 ) & ( z[j] - zbbbs <= 0 ) ) )

                #if ( q1 > 0 ) & ( q2 > 0 ) & ( q3 > 0 ) & ( q4 > 0 ) :
                #    iiInsideA[i,j]  = 1
                
        iiInside    = np.where ( iiInsideA > 0 )
        iiOutside   = np.where ( iiInsideA == 0 )

#    print nW, nH, nbbbs, limitr
#    print rdim, zdim, rcentr, rleft, zmid
#    print rmaxis, zmaxis, simag, sibry, bcentr

#   Plot output
    fig='No figure'
    if (doplot):
        N=10
        if not isinstance(doplot,bool):
            if isinstance(doplot,int):
                 N=doplot
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_aspect('equal')
            plt.contour ( r, z, psizr.T, N )
            plt.plot ( rbbbs, zbbbs, 'k', linewidth = 3 )
            if (dolimiter):
                plt.plot ( rlim, zlim, 'g', linewidth = 4 )
            plt.show ()
        else:
            ax.contour (r, z, psizr.T, N )
            ax.plot ( rbbbs, zbbbs, 'k', linewidth = 3 )
            if (dolimiter):
                ax.plot ( rlim, zlim, 'g', linewidth = 4 ) 

    #checks
    # rmaxis =/ rcentr
    eqdsk = {'nW':nW, 'nH':nH, 'nV':nV, 'nbbbs':nbbbs, 'limitr':limitr, 'rdim':rdim,
             'zdim':zdim, 'rcentr':rcentr, 'rleft':rleft, 'zmid':zmid, 
             'rmaxis':rmaxis, 'zmaxis':zmaxis, 'simag':simag, 'sibry':sibry,
             'bcentr':bcentr, 'current':current, 'fpol':fpol, 'pres':pres,
             'ffprim':ffprim, 'pprime':pprime, 'psizr':psizr, 'qpsi':qpsi, 'rbbbs':rbbbs,
             'zbbbs':zbbbs, 'rlim':rlim, 'zlim':zlim, 'r':r, 'z':z,
             'fluxGrid':fluxGrid, 'iiInside':iiInside, 'cocos':cocos, 'name':filename}

    return eqdsk,fig


def readGEQDSK2(filename='eqdsk.dat', dointerior=False, width=9, cocos=3, 
                doplot=None, dolimiter=None, ax=None, dodebug=False):
    """
    Read an eqdsk file for various cocos conventions, optionally produce a plot
    dointerior returns list of i,j pts inside LCF.
    cocos=3 is EFIT
    ax is a figure handle to plot eq on

    """
    
    import re
    import numpy as np
    import pylab as plt
    import fortranformat as ff


    f2000=ff.FortranRecordReader('a48,3i4')
    f2020=ff.FortranRecordReader('5e16.9')
    f2022=ff.FortranRecordReader('2i5')
    xdum = np.zeros(5)
    
    def readVar(fmt,line):
        return fmt.read(line)


    def readArray(fmt,shp):
        vals=[]
        if len(shp)==1: N=shp[0]
        if len(shp)==2: N=shp[0]*shp[1]
        nlines = int(N/5)
        if (N%5)!=0: nlines+=1        
        for i in range( nlines ):
            vals.extend(fmt.read(next(f)))
        return np.reshape(np.array(vals[0:N]),shp)            

        
    with open(filename, "r") as f:
        [casestr, idum, nw, nh]            =f2000.read(next(f))
        [rdim,zdim,rcentr,rleft,zmid]      =f2020.read(next(f))
        [rmaxis,zmaxis,simag,sibry,bcentr] =f2020.read(next(f))
        [current,simag,xdum,rmaxis,xdum]   =f2020.read(next(f))
        [zmaxis,xdum,sibry,xdum,xdum]      =f2020.read(next(f))
        fpol    =readArray(f2020,[nw])
        pres    =readArray(f2020,[nw])
        ffprim  =readArray(f2020,[nw])
        pprime  =readArray(f2020,[nw])
        psizr   =readArray(f2020,[nw,nh]).T # ff follows fortran indexing convention, so transpose to be consistent with usage
        qpsi    =readArray(f2020,[nw])
        #check if bb present
        [nbbbs,limitr]=f2022.read(next(f))
        RZbnd   =readArray(f2020,[nbbbs*2]) #Rbnd,Zbnd)
        RZlim   =readArray(f2020,[limitr*2]) #Rlim,Zlim)

    if dodebug: 
        print("Data string header:", rdim,zdim,rcentr,rleft,zmid )
        print("Dimensions:", nw, nh, nbbbs, limitr, rdim, zdim )

    rbbbs,zbbbs  = RZbnd[::2],RZbnd[1::2]
    rlim,zlim    = RZlim[::2],RZlim[1::2]

    if dodebug: print('one D arrays: ', fpol[-1],pres[-1], ffprim[-1], pprime[-1], qpsi[-1] )

    rStep   = rdim / ( nw - 1 )
    zStep   = zdim / ( nh - 1 )
    fStep   = -( simag - sibry ) / ( nw - 1 )

    r   = np.arange ( nw ) * rStep + rleft
    z   = np.arange ( nh ) * zStep + zmid - zdim / 2.0

    fluxGrid    = np.arange ( nw ) * fStep + simag

#   Find indices of points inside and outside
#   the rbbbs/zbbbs boundary.
    import matplotlib.path as mplPath
    import numpy as np
    lcf=mplPath.Path( np.column_stack( (rbbbs,zbbbs) ) )
    iiInsideA   = np.zeros ( psizr.shape )
    iiInside = -1
    iiOutside = -1
    if (dointerior):
        for i in np.arange ( nw ) :
            for j in np.arange ( nh ) :
                if lcf.contains_point( (r[i],z[i]) ):
                    iiInsideA[i,j] = 1
                
        iiInside    = np.where ( iiInsideA > 0 )
        iiOutside   = np.where ( iiInsideA == 0 )

#   Plot output
    fig='No figure'
    if (doplot):
        N=10
        if not isinstance(doplot,bool):
            if isinstance(doplot,int):
                 N=doplot
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_aspect('equal')
            plt.contour ( r, z, psizr.T, N )
            plt.plot ( rbbbs, zbbbs, 'k', linewidth = 3 )
            if (dolimiter):
                plt.plot ( rlim, zlim, 'g', linewidth = 4 )
            plt.show ()
        else:
            ax.contour (r, z, psizr.T, N )
            ax.plot ( rbbbs, zbbbs, 'k', linewidth = 3 )
            if (dolimiter):
                ax.plot ( rlim, zlim, 'g', linewidth = 4 ) 


    eqdsk = {'nW':nw, 'nH':nh, 'nbbbs':nbbbs, 'limitr':limitr, 'rdim':rdim,
             'zdim':zdim, 'rcentr':rcentr, 'rleft':rleft, 'zmid':zmid, 
             'rmaxis':rmaxis, 'zmaxis':zmaxis, 'simag':simag, 'sibry':sibry,
             'bcentr':bcentr, 'current':current, 'fpol':fpol, 'pres':pres,
             'ffprim':ffprim, 'pprime':pprime, 'psizr':psizr, 'qpsi':qpsi, 'rbbbs':rbbbs,
             'zbbbs':zbbbs, 'rlim':rlim, 'zlim':zlim, 'r':r, 'z':z, 'psirz':psizr.T,
             'fluxGrid':fluxGrid, 'iiInside':iiInside, 'cocos':cocos, 'name':filename}

    return eqdsk,fig


def getModB(eq):
    """
    Calculate the magnitude of the magnetic field on the RZ mesh.


        |B| = \\sqrt(Fpol^2+(d\\Psi/dZ)^2+(d\\Psi/dR)^2)/R

    where Fpol== R*Bphi , Bpol = |grad Psi|/R
    """
    import numpy as np
    from scipy import interpolate

    #poloidal component. for cocos=[3] or 1/11 (R,phi,Z)
    fluxfactor=1.0 ; sbp = +1
    if eq['cocos']>=11   : fluxfactor=2.*np.pi
    if eq['cocos']%10==3 : sbp=-1
    
    R=eq.get('r')
    Z=eq.get('z')
    Rv,Zv=np.meshgrid(R,Z,indexing='ij') #these are R and Z on RZ mesh, first index for Z , default indexing
    psiRZ=eq.get('psizr').T 
    psiZR=eq.get('psizr')   
    spline_psi = interpolate.RectBivariateSpline(R,Z,psiRZ.T,bbox=[np.min(R),np.max(R),np.min(Z),np.max(Z)],kx=5,ky=5)
    psi_int_r=spline_psi.ev(Rv,Zv,dx=1)/fluxfactor
    psi_int_z=spline_psi.ev(Rv,Zv,dy=1)/fluxfactor
    grad_psi=np.sqrt(psi_int_z**2+psi_int_r**2)
    
    #toroidal component
    #get Fpol and interpolate to RZ mesh to get fpolRZ
    fpol=eq.get('fpol')
    psi=eq.get('fluxGrid') #mesh for fpol
    #fpolRZ=fpol(psiRZ)
    # scipy 0.18    spline_fpol=interpolate.CubicSpline(psi,fpol,bc_type='natural')
    #fill value is B0*R0
    spline_fpol=interpolate.interp1d(psi,fpol,bounds_error=False,fill_value=fpol[-1],kind='cubic')
    fpolRZ=[] #np.zeros(psiRZ.shape)
    for psirow in psiZR:
        fpolRZ.append( spline_fpol(psirow) )
    fpolRZ=np.array(fpolRZ) #Fpol numpy array on RZ mesh

    modgradpsi=np.sqrt(grad_psi**2+fpolRZ**2)
    modB=modgradpsi/Rv
    if R[0]==0.0: #If origin is included in domain, be careful with |B| on axis.
        modB[:,0]=(np.diff(modgradpsi,axis=1)/(R[1]-R[0]))[:,0]
    BV=( +sbp*psi_int_z/Rv, fpolRZ/Rv, -sbp*psi_int_r/Rv) #R,phi,Z for cocos1/11 and 3/13
    #Add components
    return modB,grad_psi,fpolRZ,Rv,Zv,BV


def getLCF(eq):
    #find which contour in LCF, same as rbbbs? 

    import matplotlib.path as mplPath
    import numpy as np
    import pylab as p

    R=eq.get('r')
    Z=eq.get('z')
    psiRZ=np.transpose(eq.get('psizr'))    
    CSlcf=p.contour(R,Z,psiRZ,levels=[eq['sibry']-.01])
    cntr=(eq['rmaxis'],eq['zmid'])
    lcf=(0,0)
    for p in CSlcf.collections[0].get_paths():
        v = p.vertices
        x = v[:,0]
        y = v[:,1]
        bbPath = mplPath.Path(np.column_stack( (x,y)))
        if bbPath.contains_point(cntr):
            lcf=(x,y)
            return lcf

        
def plotEQDSK(eq):
    import pylab as plt
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(10, 6))
    fig.suptitle( 'EQDSK content for '+eq['name'] )
    
    modB,grad_psi,fpolRZ,Rv,Zv,BV=getModB(eq)
    R=eq.get('r')
    Z=eq.get('z')
    nr2=int(eq['nW']/2)
    nz2=int(eq['nH']/2)


    ax1.grid()
    ax1.set_title("Magnetic field components")
    ax1.plot(R,modB [:,nz2-1],'purple',label='B')
    ax1.plot(R,BV[1][:,nz2-1],'black',label='Btor')
    ax1.plot(R,BV[2][:,nz2-1],'orange',label='BZ')
    ax1.plot(R,BV[0][:,nz2-1],'red',label='BR')
    ax1.legend()

    ax2.set_title('Flux surfaces')
    ax2.contour (eq['r'], eq['z'], eq['psizr'].T, 40 )
    ax2.plot ( eq['rbbbs'], eq['zbbbs'], 'k', linewidth = 3 )
    ax2.plot ( eq['rlim'], eq['zlim'], 'g', linewidth = 4 )
    ax2.set_aspect('equal')
    
    ax3.set_title('Profiles')
    ax3.plot(eq['fluxGrid'], eq['qpsi'],                   label='q')
    ax3.plot(eq['fluxGrid'], eq['fpol']/ eq['fpol'][0],    label='F/F(0)')
    ax3.plot(eq['fluxGrid'], eq['pres']/eq['pres'][0],     label='p/p(0)')
    ax3.plot(eq['fluxGrid'], eq['ffprim']/eq['ffprim'][0], label="FF' norm")
    ax3.plot(eq['fluxGrid'], eq['pprime']/eq['pprime'][0], label="p' norm")    
    ax3.legend()
    
    ax4.text(0.5,0.9,'Values from EQDSK header.',ha='center')
    hcol=0
    for i, (key, value) in enumerate(eq.items()):
        if i>15: break
        if i>8: hcol=1
        ax4.text(0.1+hcol*0.4, 0.9 - (i + 1) * 0.09+hcol*8*0.09,
                 f'{key}: {value:5.2f}', ha='left', va='center', fontsize=12)


def writeEQDSK(eq,fname):
    """Write out the equilibrium in G-EQDSK format.
    Code courtesy of eq_JET.py from OMFIT"""

    import fortranformat as ff
    import numpy as np

    # Open file for writing
    f = open(fname, 'w')

    nr = eq['nW']
    nz = eq['nH']    

    # Get eq at this timeslice
    rdim    = eq['rdim']
    zdim    = eq['zdim']
    rcentr  = eq['rcentr']
    rleft   = eq['rleft']
    zmid    = eq['zmid']
    rmaxis  = eq['rmaxis']
    zmaxis  = eq['zmaxis']
    simag   = eq['simag']
    sibry   = eq['sibry']
    bcentr  = eq['bcentr']
    current = eq['current']
    xdum    = 0.0

    def GetSlice( data, N, ti ):
        return data[ ti * N : ( ti + 1 ) * N ]

    # FPOL eq
    fpol =  eq['fpol']

    # Pressure eq
    pressure = eq['pres']

    # FFPRIM eq
    ffprim = eq['ffprim']

    # PPRIME eq
    pprime = eq['pprime']

    # PSI eq
    psi = np.transpose(eq['psizr'])

    # Q eq
    q = eq['qpsi']
    
    # Plasma Boundary
    Rbnd = eq['rbbbs']
    Zbnd = eq['zbbbs']
    n_bnd = eq['nbbbs'] 

    # Limiter eq
    Rlim = eq['rlim']
    Zlim = eq['zlim']
    limitr = len(Rlim)

    # Write Eqdsk from -----------------------------------

    f2020=ff.FortranRecordWriter('5e16.9')
    f2022=ff.FortranRecordWriter('2i5')
        
    def writeVar(handle,varList):
        f.write(handle.write(varList))
        f.write("\n")

    def writeOrderedPairs(handle,var1,var2):
        longArrayOfPairs=[]
        for pv,_ in enumerate(var1):
            longArrayOfPairs.append(var1[pv])
            longArrayOfPairs.append(var2[pv])

        writeVar(handle,longArrayOfPairs)
        
    A52 = 'plasma ep_v2.0_:_01:01:17'.ljust(48)
    f.write(A52[0:48])
    writeVar(ff.FortranRecordWriter('3i4'), [0,nr,nz] )
    writeVar(f2020,[rdim,zdim,rcentr,rleft,zmid])
    writeVar(f2020,[rmaxis,zmaxis,simag,sibry,bcentr])
    writeVar(f2020,[current,simag,xdum,rmaxis,xdum]) 
    writeVar(f2020,[zmaxis,xdum,sibry,xdum,xdum])  
    writeVar(f2020,fpol)
    writeVar(f2020,pressure)
    writeVar(f2020,ffprim)
    writeVar(f2020,pprime)
    writeVar(f2020,psi.flatten())
    writeVar(f2020,q)
    writeVar(f2022,[n_bnd,limitr])
    writeOrderedPairs(f2020,Rbnd,Zbnd)
    writeOrderedPairs(f2020,Rlim,Zlim)
    
    f.close()


def resize(nx,eq):
    import copy
    neweq=copy.deepcopy(eq)
    neweq['nW']=nx
    neweq['nH']=nx

    def resize1D(x,newx,profile):
      import numpy as np
      from scipy import interpolate
      f = interpolate.interp1d(x, profile)
      return f(newx)
      

    nW=nx
    nH=nx
    rdim=eq['rdim']
    zdim=eq['zdim']
    simag=eq['simag']
    sibry=eq['sibry']
    rStep   = rdim / ( nW - 1 )
    zStep   = zdim / ( nH - 1 )
    fStep   = -( simag - sibry ) / ( nW - 1 )
    r   = np.arange ( nW ) * rStep + rleft
    z   = np.arange ( nH ) * zStep + zmid - zdim / 2.0
    fluxGrid    = np.arange ( nW ) * fStep + simag

    neweq['r']=r
    neweq['z']=z
    neweq['fluxGrid']=fluxGrid

    neweq['fpol']=resize1D(eq['fluxGrid'],fluxGrid,eq['fpol'])
    neweq['pres']=resize1D(eq['fluxGrid'],fluxGrid,eq['pres'])
    neweq['ffprim']=resize1D(eq['fluxGrid'],fluxGrid,eq['ffprim'])
    neweq['pprime']=resize1D(eq['fluxGrid'],fluxGrid,eq['pprime'])

    return neweq


def rescaleB(eq,filename,s=1.,sR=1.):
    """
    Rescale by MHD scalings.
    Psi scaling: Psi=s*Psi, P=s^2 P, g=s*g, beta stays the same
    Pressure scaling: P=P+c, Psi unchanged
    Toroidal field scaling: g^2=g^2+c, Psi unchanged
    """
    import copy

    R0=eq['rmaxis']
    f=newR/R0
    f=sR

    neweq=copy.deepcopy(eq)
    neweq['psizr']=eq['psizr']*f
    neweq['fpol']=eq['fpol']*f

    neweq['bcentr']=eq['bcentr']*f
    neweq['ffprim']=eq['ffprim']*f*f
    neweq['simag']=eq['simag']*f    
    neweq['sibry']=eq['sibry']*f 
    neweq['current']=eq['current']*f 
    neweq['fluxgGrid']=eq['fluxGrid']*f

    writeEQDSK(neweq,filename)


