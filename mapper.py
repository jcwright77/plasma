import scipy.integrate as integrate
from scipy.interpolate import griddata
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator
import numpy as np
#ensure periodicity with fft
import scipy.fftpack as sft
from scipy import integrate
import scipy.interpolate

from plasma import  equilibrium_process as eqdsk

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.path import Path
from packaging.version import Version

mu0=4.*np.pi*1.e-7
def calcarea(x,y):
    #x=vs[:,0]
    #y=vs[:,1]
    return 0.5*np.sum(y[:-1]*np.diff(x) - x[:-1]*np.diff(y))

def mapper(eqobj,jac='eqarc'):
  """
    mapper calculates a r,theta cooridinate system within the last closed
    flux surface
    eqfile: filename of geqdesk file or dictionary of values from geqdsk file
    jac: 'straight' or 'eqarc' 

    returnx Xmap,Zmap on r,theta grid
  """
      
  dodebug=False #True
  sepfrac=0.995  
  #print("type",str(type(eqobj)))
  if isinstance(eqobj,str):
    eq=eqdsk.readGEQDSK(eqobj)[0]
  elif isinstance(eqobj,dict):
    eq=eqobj
  else:
    return "Unrecognized equilibrium object, must be filename or dictionary"

  R=eq.get('r')
  Z=eq.get('z')-eq['zmaxis']  #axis needs to be at z=0 for mapping
  B,grad_psi,fRZ,Rv,Zv,Bv=eqdsk.getModB(eq)
  psi=eq.get('psizr').T
  curtor=[]
  area=[]
  mapzmaxis=float(0.0)
  ifrhopol=True #False #use root psipol mesh instead of psipol (eg for torlh)


  def find_cut(x,y, rmaxis, zmaxis):
    #adapted from S. Shiraiwai to find crossing going counter clockwise
    for k in range(len(y)-1):
      km = k-1
      if y[km] < zmaxis and y[k] > zmaxis:
        return k
    return -1


  nsample=1200  
  r200=np.linspace(min(R),max(R),1200)
  z200=np.linspace(min(Z),max(Z),1200)
  RR,ZZ=np.mgrid [min(R):max(R):1200j, min(Z):max(Z):1200j ]

  spline_psi = scipy.interpolate.RectBivariateSpline(R,Z,psi.T,bbox=[np.min(R),
                                      np.max(R),np.min(Z),np.max(Z)],kx=5,ky=5)
  psi_int=spline_psi.ev(RR,ZZ)
  psi_int_r=spline_psi.ev(RR,ZZ,dx=1)
  psi_int_z=spline_psi.ev(RR,ZZ,dy=1)
  grad_psi=np.sqrt(psi_int_z**2+psi_int_r**2)

  #generated mapped mesh size:
  npsi=80
  ntheta=128

  #we will have eq_theta at each filtered_cx,cy
  #like polar contour needed to be converted to regular mesh

  #get psi mesh for surfaces for theta within LCF
  #LCS is at psi=0
  #drop first point, magnetic axis. We add this one manually since it cannot
  #be contoured.

  #psimesh=eq['fluxGrid'][1:] #poloidal flux grid, created by readGEQDSK from eqdsk but not in eqdsk file
  #resize psi to the number of desired psi levels, psimesh is already uniform
  #initial psimesh is [-psimin,0].
  #the following is only necessary if psimesh is not uniform which it should be for an eqdsk file.

  #nidx=100
  sgnflux=np.sign(  eq['simag']+eq['sibry']  )
  if ifrhopol:
    rhopol = np.sqrt(np.linspace( np.abs(eq['simag']),np.abs(eq['sibry'])*sepfrac,npsi)*sgnflux)
    fity = rhopol**2*sgnflux #reference psipol consistent with uniform rhopol
    rhopol = np.linspace(0,1,npsi)
    psimesh=fity
    eq['rhopolmap']=rhopol  #sqrt norm rho pol for map size npsi, linear spaced
  else:
    psimesh=np.linspace(eq['simag'],eq['sibry']*sepfrac,npsi)
    eq['rhopolmap']=np.sqrt(np.linspace(0.,sepfrac,npsi))
    rhopol = np.linspace(0,1,npsi)

  rmaxis = eq['rmaxis']
  zmaxis = mapzmaxis #eq['zmaxis']
  eq['psipolmap'] = psimesh
  if dodebug: print('psimesh',psimesh)
  if dodebug: print('psiaxis',eq.get('simag'))
  c_pprime  = np.interp( psimesh, eq['fluxGrid'], eq['pprime'] )
  c_ffprime = np.interp( psimesh, eq['fluxGrid'], eq['ffprim'] )
  qmap      = np.interp( psimesh, eq['fluxGrid'], eq['qpsi']   )
  #calculate toroidal flux and rhotor
  psitormap=integrate.cumulative_trapezoid(qmap,psimesh,initial=0.0)
  rhotor=(psitormap-psitormap[0])/(psitormap[-1]-psitormap[0])
  eq['rhotormap']=rhotor
    
#Extract contours and values for flux coordinate system.
#contours go counter-clockwise, which we want
#contours don't necessarily start at y=0., so rebase
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_aspect('equal')
  fig.set_figheight(6)
  psi_cs=plt.contour(RR,ZZ,psi_int,levels=psimesh)

  psixy=[]
  if dodebug: print('matplotlib version', str(Version(matplotlib.__version__)))
#Get Psi contours, careful to exclude field coils
  if Version(matplotlib.__version__)  < Version('3.8.0'):
    for p in psi_cs.collections:
      for pp in p.get_paths():
        v = pp.vertices
        x = v[:,0]
        y = v[:,1]
        if np.abs(np.average(y))<0.10*eq['rmaxis'] and np.abs(np.average(x))<0.10*eq['rmaxis']: #only keep core plasma contours
          psixy.append( (x,y) )
  else:
    for i,crvs in enumerate(psi_cs.allsegs):
      knds=psi_cs.allkinds[i]
      for j,crv in enumerate(crvs):
        x,y=zip(*crv)
        knd=knds[j]
        hasaxis=Path(crv,knd).contains_point( (eq['rmaxis'],eq['zmaxis'])  )      
        if hasaxis:
          psixy.append( (x,y) )


#Define uniform theta mesh
  uni_theta=np.linspace(0,2.0*np.pi,ntheta,endpoint=False)

#Set up X(psi,theta) Y(psi,theta) and initialize with magnetic axis point
  eq_x=[]
  eq_y=[]

#  print('grad psi, B shapes', grad_psi.shape, B.shape)
  #points = np.array( (RR.flatten(), ZZ.flatten()) ).T
  gpsivalues = grad_psi.flatten()
  spline_psi = scipy.interpolate.RectBivariateSpline(r200,z200,grad_psi)
  #Bpoints = np.array( (RR.flatten(), ZZ.flatten()) ).T
  Bvalues   = B.flatten()
  rB = np.linspace(min(R),max(R),B.shape[0])
  zB = np.linspace(min(R),max(R),B.shape[1])
  spline_B = scipy.interpolate.RectBivariateSpline(rB,zB,B)

  minmodes = 5
  maxmodes = 12
  for c_idx,(cx,cy) in enumerate(psixy):
      #for each surface, low pass filter to central 8+ DC Fourier modes
      #remove last element for fft since it is equal to the first element
      #size of (cx,cy) is  variable
      #print('cx',c_idx, len(cx),len(cy),len(psixy),len(psixy[0]),psimesh[c_idx]) #just to see progress
      area1=area2
      area2=calcarea(cx,cy)
      #shift to midplane as first element
      idx=find_cut(cx,cy,rmaxis,zmaxis)
      if (idx>=0):
          #now project to centers
          #this assumes that values stradle midplane which seems to be the
          #case for python contour, but should be made more robust.
          cx=0.5*(np.roll(cx,-idx)+np.roll(cx,-idx+1))
          cy=0.5*(np.roll(cy,-idx)+np.roll(cy,-idx+1))

      #filter out high freq noise, esp needed near axis
      nmodes=max( minmodes,int(len(cx)/float(maxmodes) ))
      if dodebug: print('c_index,nmodes',c_idx,nmodes)

      fftx=sft.fft(cx)
      fftx[int(nmodes/2)+1:-int(nmodes/2)]=0
      filtered_cx=sft.ifft(fftx).real

      ffty=sft.fft(cy)
      ffty[int(nmodes/2)+1:-int(nmodes/2)]=0

      #restore value for idx=0 for Y.
      ffttotal=np.sum(ffty) #want to restore to total before filter.
      ffty[int(nmodes/2)+1]=-ffttotal/2.
      ffty[-int(nmodes/2)-1]=-ffttotal/2.

      filtered_cy=sft.ifft(ffty).real

      #interpolate |grad psi| and B onto this surface
      #this *significantly* slows down this routine. 

      # filtered_cx=cx ; filtered_cy=cy
      if jac=="straight":
          c_B  = spline_B.ev(filtered_cx,filtered_cy)
         # c_B  = griddata( Bpoints, Bvalues, (filtered_cx,filtered_cy), method='cubic' )


      #c_gradpsi = griddata( points, gpsivalues, (filtered_cx,filtered_cy), method='cubic' )
      c_gradpsi  = spline_psi.ev(filtered_cx,filtered_cy)
      
      #derivative from fft needs factor of 2pi
      df_dx=sft.diff(filtered_cx)*np.pi*2.0/len(filtered_cx)
      df_dy=sft.diff(filtered_cy)*np.pi*2.0/len(filtered_cy)
      dl=np.sqrt(df_dx**2+df_dy**2) #These two steps could be done with FFT too.

      #flux surface integrals here for later use, eg current profile
      c_curtor = -integrate.simpson(
          dl/c_gradpsi* filtered_cx*( c_pprime[c_idx] +
                                      c_ffprime[c_idx]/filtered_cx**2/mu0 ) )


      #From G-S J_phi force balance equation
      c_area = integrate.simpson(
          dl/c_gradpsi
      ) #this area is centered on the cell darea/dpsi
      area.append( c_area )
      curtor.append (c_curtor/c_area)  #make this d<Jphi>/dpsi

      #if straight field line, multiply dl by 1/R*|gradpsi|
      if jac=='straight':
        dtheta=dl/np.abs(c_gradpsi*filtered_cx)
      else: #jac='eqarc'
        dtheta=dl
      L=integrate.cumulative_trapezoid(dtheta,initial=0)/len(dtheta)

      #now put each on same theta mesh, Jacobian selection
      this_theta=L/np.max(L)*2.0*np.pi
      t_map=scipy.interpolate.interp1d(this_theta,filtered_cx,kind='cubic')
      eq_x.append(t_map(uni_theta))
      t_map=scipy.interpolate.interp1d(this_theta,filtered_cy,kind='cubic')
      eq_y.append(t_map(uni_theta))

  #add origin term
  #area.insert(0,0.)          #area of origin is 0.
  #curtor.insert(0,curtor[0]) #current density maximum at origin
  #print('curtor size', len(curtor), len(area) )
  curtor=np.array(curtor)
  area=np.array(area)
  eq['darea']=area
  eq['Jtor']=curtor

  #print('curtor size', len(curtor), len(area) )
  #center_psimap=(( eq['psipolmap']+np.roll(eq['psipolmap'],1)  )/2)
  #Ipsi=scipy.integrate.cumulative_trapezoid( eq['Jtor'],eq['darea'], initial=0)
  Ipsi=scipy.integrate.cumulative_trapezoid(eq['Jtor']*eq['darea']*np.diff(eq['psipolmap']),initial=0) #),3.14159*0.01)

  print("Current accuracy eqdsk,mapper", eq['current'],Ipsi[-1]," rescaling")
  Ipsi_mod=Ipsi*eq['current']/Ipsi[-1] 
#  eq['Jtor']=eq['Jtor']*eq['current']/Ipsi[-1] 
  #Ipsi = scipy.integrate.cumulative_trapezoid( curtor, center_psimap, initial=0)
  eq['Ipsi']=Ipsi_mod
  
  #center_area=(( eq['darea_dpsi']+np.roll(eq['darea_dpsi'],1)  )/2)[1:]
  #totarea=scipy.integrate.cumulative_trapezoid( eq['darea_dpsi'][1:], center_psimap)
  
  Xmap=np.double(np.vstack(eq_x))
  Zmap=np.double(np.vstack(eq_y))
  #add origin
  NXmap=np.zeros([npsi,ntheta])
  NZmap=np.zeros([npsi,ntheta])
  #print('Xmap',Xmap.shape,NXmap.shape,len(eq_x))
  zax,rax=np.average((Zmap[0,:]-mapzmaxis)),np.average((Xmap[0,:]))
  NZmap[0,:]= 0.0
  #NZmap= NZmap - zax
  NXmap[0,:]=rax #eq['zmaxis']
  NXmap[1:,:]=Xmap
  NZmap[1:,:]=Zmap

  print("Error or shift of vertical axis from zero is",
        np.average((Zmap[0,:]-mapzmaxis)),np.average((Xmap[0,:])),
        eq['zmaxis'],eq['rmaxis']
       )
  print(Zmap[0,:],mapzmaxis)
  print("Error or shift of vertical axis from zero is",
        np.average((Zmap[1,:]-mapzmaxis)),np.average((Xmap[1,:])),
        eq['zmaxis'],eq['rmaxis']
       )
  del(Xmap)
  del(Zmap)
  Xmap=NXmap
  Zmap=NZmap
  #Zmap[0,0]=np.double(0.0)

  print('shapes of mapped arrays: ',Xmap.shape,Zmap.shape)
  #save this in pickle file
  eq['xmap']=Xmap
  eq['zmap']=Zmap
  eq['jac']=jac
  return eq


def plot_equilibrium(eq):

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_aspect('equal')
  fig.set_figheight(4)

  Xmap=eq['xmap'] ; Zmap=eq['zmap']
  maxpsi=0.995
  maxpsiind=int(maxpsi*Xmap.shape[0])
  #Theta lines
  for i in np.arange(0,len(Xmap[0,:]),5):
    ax.plot(Xmap[:maxpsiind,i],Zmap[:maxpsiind,i])

  #Psi surface
  for i in np.arange(0,len(Xmap[:maxpsiind,0]),5):
    ax.plot(Xmap[i,:],Zmap[i,:])
  ax.set_title('Surfaces of constant theta and psi (every 10th)');
