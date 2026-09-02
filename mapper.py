import numpy as np

import scipy.integrate as integrate
from scipy.interpolate import griddata
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator
import scipy.fftpack as sft
from scipy import integrate
import scipy.interpolate

from plasma import  equilibrium_process as eqdsk

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.path import Path
import matplotlib.patches as mpatches
from packaging.version import Version

import copy

mu0=4.*np.pi*1.e-7
def calcarea(x,y):
    #x=vs[:,0]
    #y=vs[:,1]
    return 0.5*np.sum(y[:-1]*np.diff(x) - x[:-1]*np.diff(y))


def min_in_polygon(X, Y, Z, px, py, findmax=False, doplot=False):
    """
    Find the minimum value of a 2D field inside a polygon.
 
    Parameters
    ----------
    X : array_like, shape (nx,)
        1D x-coordinates of the rectangular mesh.
    Y : array_like, shape (ny,)
        1D y-coordinates of the rectangular mesh.
    Z : array_like, shape (ny, nx)
        2D field values; Z[iy, ix] corresponds to point (X[ix], Y[iy]).
    px : array_like, shape (n,)
        x-coordinates of the polygon vertices.
    py : array_like, shape (n,)
        y-coordinates of the polygon vertices.
 
    Returns
    -------
    max_val : float
        Maximum value of Z inside the polygon.
    ix_max : int
        Column index (into X) of the maximum.
    iy_max : int
        Row index (into Y) of the maximum.
    x_max : float
        x-coordinate of the maximum.
    y_max : float
        y-coordinate of the maximum.
 
    Raises
    ------
    ValueError
        If no mesh points fall inside the polygon.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    Z = np.asarray(Z, dtype=float)
    px = np.asarray(px, dtype=float)
    py = np.asarray(py, dtype=float)
 
    if Z.shape != (len(Y), len(X)):
        raise ValueError(
            f"Z shape {Z.shape} does not match (len(Y), len(X)) = ({len(Y)}, {len(X)})"
        )
 
    # --- Build the polygon path (auto-close) ----------------------------
    # matplotlib Path expects the polygon to be explicitly closed
    verts = np.column_stack([px, py])
    if not np.array_equal(verts[0], verts[-1]):
        verts = np.vstack([verts, verts[0]])  # close the ring
    poly = Path(verts)
 
    # --- Create a grid of all mesh points --------------------------------
    # Meshgrid: XX[iy, ix] = X[ix], YY[iy, ix] = Y[iy]
    XX, YY = np.meshgrid(X, Y)              # both shape (ny, nx)
    points = np.column_stack([XX.ravel(), YY.ravel()])  # shape (ny*nx, 2)
 
    # --- Mask: True where the point is inside the polygon ----------------
    inside = poly.contains_points(points)   # shape (ny*nx,)
    inside_2d = inside.reshape(Z.shape)     # shape (ny, nx)
 
    if not inside_2d.any():
        raise ValueError("No mesh points found inside the polygon.")
 
    # --- Masked array and argmax -----------------------------------------
    if findmax:
        Z_masked = np.where(inside_2d, Z, -np.inf)
        flat_idx  = np.argmax(Z_masked)
    else:
        Z_masked = np.where(inside_2d, Z, +np.inf)
        flat_idx  = np.argmin(Z_masked)
    iy_max, ix_max = np.unravel_index(flat_idx, Z.shape)
 
    max_val = Z[iy_max, ix_max]
    x_max   = X[ix_max]
    y_max   = Y[iy_max]

    if doplot:
        fig, ax = plt.subplots(figsize=(7, 6))
        c = ax.contour(X, Y, Z, levels=40, cmap="viridis")
        fig.colorbar(c, ax=ax, label="Z value")
        ax.set_aspect('equal')
 
        poly_closed = np.append(px, px[0]), np.append(py, py[0])
        ax.plot(*poly_closed, "k-", lw=2, label="Polygon")
        ax.plot(x_max, y_max, "r*", markersize=18, label=f"Max = {max_val:.3f}")
 
        ax.set_title("Minimum in 2-D field inside polygon")
        ax.legend(loc="upper left")
        plt.tight_layout()
        
    return max_val, ix_max, iy_max, x_max, y_max
 
 
def mapper(eqobj,jac='eqarc',maxmom=12, npsi=40, ntheta=128, nsample=600,
           sepfrac=0.98,dodebug=False,doplot=False,ifrhopol=True):
  """
    mapper calculates a r,theta cooridinate system within the last closed
    flux surface
    eqfile: filename of geqdesk file or dictionary of values from geqdsk file
    jac: 'straight' or 'eqarc' 

    returns Xmap,Zmap on r,theta grid
  """
      
  if isinstance(eqobj,str):
    eq=eqdsk.readGEQDSK(eqobj)[0]
  elif isinstance(eqobj,dict):
    eq=copy.deepcopy(eqobj)
  else:
    return "Unrecognized equilibrium object, must be filename or dictionary"

  R=eq.get('r')
  Z=eq.get('z')#-eq['zmaxis']  #axis needs to be at z=0 for mapping
  B,grad_psi,fRZ,Rv,Zv,Bv=eqdsk.getModB(eq)
  psi=eq.get('psizr').T
  curtor=[]
  area=[]
  mapzmaxis=float(0.0)


  def find_cut(x,y, rm, zm):
    #adapted from S. Shiraiwa to find crossing going counter clockwise
    for k in range(len(y)-1):
      km = k-1
      if y[km] < zm and y[k] > zm:
        return k
    return -1

  r200=np.linspace(min(R),max(R),nsample)
  z200=np.linspace(min(Z),max(Z),nsample)
  RR,ZZ=np.mgrid [min(R):max(R):np.complex64(0,nsample), min(Z):max(Z):np.complex64(0,nsample) ]

  spline_psi = scipy.interpolate.RectBivariateSpline(R,Z,psi.T,bbox=[np.min(R),
                                      np.max(R),np.min(Z),np.max(Z)],kx=5,ky=5)
  psi_int=spline_psi.ev(RR,ZZ)
  psi_int_r=spline_psi.ev(RR,ZZ,dx=1)
  psi_int_z=spline_psi.ev(RR,ZZ,dy=1)
  grad_psi=np.sqrt(psi_int_z**2+psi_int_r**2)


#Define uniform theta mesh
  uni_theta=np.linspace(0,2.0*np.pi,ntheta,endpoint=False)

#Set up X(psi,theta) Y(psi,theta) and initialize with magnetic axis point
  eq_x=[]
  eq_y=[]

  rmaxis = eq['rmaxis']
  zmaxis = eq['zmaxis']
  
  spline_gpsi = scipy.interpolate.RectBivariateSpline(r200,z200,grad_psi)
  
  #check axis position
  max_val, ix_max, iy_max, x_max, y_max = min_in_polygon(r200, z200, (psi_int.T),
                                                         eq['rlim'], eq['zlim'],doplot=doplot)
  if dodebug: print('maxind2',  max_val, ix_max, iy_max, x_max, y_max,rmaxis,zmaxis)
  eq['rmaxis'] = x_max
#  eq['zmaxis'] = y_max
  rmaxis=x_max #; zmaxis=y_max
  
  rB = np.linspace(min(R),max(R),B.shape[0])
  zB = np.linspace(min(R),max(R),B.shape[1])
  spline_B = scipy.interpolate.RectBivariateSpline(rB,zB,B)

  
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

  sgnflux=np.sign( eq['simag']+eq['sibry']  )
  if eq['simag']>eq['sibry'] :
      print('Error, mapper requires increase poloidal flux, please convert to cocos%10=1,2,5,6 first')
      exit

  dpsi=(eq['sibry']-eq['simag'])/(npsi-2.)
  simax=(eq['sibry']-eq['simag'])*sepfrac+eq['simag']
  simin=eq['simag']+dpsi
  if ifrhopol:  #su btract 1 from npsi to add origin later but not try to contour it
    sgnpsi = np.sign(np.linspace( simin,simax,npsi-1) )
    rhopol = np.linspace( np.sqrt(np.abs(simin)),np.sqrt(np.abs(simax) ),npsi-1)
    fity = rhopol**2*sgnpsi #values of flux space uniformly approx in space
    rhopol = np.linspace(0,1,npsi) #rhopol is just 0,1 mesh uniform
    psimesh=fity
    eq['rhopolmap']=rhopol  #sqrt norm rho pol for map size npsi, linear spaced
  else:
    psimesh=np.linspace(simin,simax,npsi)
    eq['rhopolmap']=np.sqrt(np.linspace(0.,sepfrac,npsi-1))
    rhopol = np.linspace(0,1,npsi)


  eq['psipolmap'] = psimesh
  if dodebug: print('psimesh',psimesh)
  if dodebug: print('rho sizes',len(psimesh),len(rhopol))
  if dodebug: print('psiaxis',eq.get('simag'),eq.get('sibry'),sgnflux)
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
        #only keep core plasma contours
        if np.abs(np.average(y))<0.05*rmaxis and np.abs(np.average(x))<0.10*rmaxis:
          psixy.append( (x,y) )
  else:
    if dodebug: print('#segs', len(psi_cs.allsegs) ,len(psimesh) )
    for i,crvs in enumerate(psi_cs.allsegs):
      knds=psi_cs.allkinds[i]
      for j,crv in enumerate(crvs):
        if dodebug: print('crv',i,j,len(psi_cs.allsegs),crv.shape)
        x,y=zip(*crv)
        knd=knds[j]
        hasaxis=Path(crv,knd).contains_point( (rmaxis,zmaxis)  ) 
        if hasaxis: #this includes the axis point in contour
          psixy.append( (x,y) )


  minmodes = 5
  maxmodes = maxmom
  area1=0.
  area2=0.  
  for c_idx,(cx,cy) in enumerate(psixy):
      #for each surface, low pass filter to central 8+ DC Fourier modes
      #remove last element for fft since it is equal to the first element
      #size of (cx,cy) is  variable

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
      nmodes=max( minmodes,int(len(cx)/float(maxmodes)/8. ))
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

      filtered_cx=cx ; filtered_cy=cy
      if jac=="straight":
          c_B  = spline_B.ev(filtered_cx,filtered_cy)

      c_gradpsi  = spline_gpsi.ev(filtered_cx,filtered_cy)
      
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
      ) #this area is centered on the cell and is darea/dpsi
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

  curtor=np.array(curtor)
  area=np.array(area)
  eq['darea']=area
  eq['Jtor']=curtor
  midpsimap= (eq['psipolmap'][1:] + eq['psipolmap'][:-1])/2
  
  if dodebug: print('curtor size', len(curtor), len(area), len(midpsimap),
                    len(np.diff(eq['psipolmap'])), len(eq['psipolmap']),len(psixy) )
  #center_psimap=(( eq['psipolmap']+np.roll(eq['psipolmap'],1)  )/2)
  #Ipsi=scipy.integrate.cumulative_trapezoid( eq['Jtor'],eq['darea'], initial=0)
#  Ipsi=scipy.integrate.cumulative_trapezoid(eq['Jtor']*eq['darea']*np.diff(eq['psipolmap']),initial=0) #),3.14159*0.01)
  #add origin pt

  Ipsi=scipy.integrate.cumulative_trapezoid( eq['Jtor']*eq['darea'], psimesh,initial=0)

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
  if dodebug: print('Xmap shape',Xmap.shape,NXmap.shape)
  
  NZmap[0,:]= zmaxis #could be vertically shifted
  NXmap[0,:]= rmaxis
  NXmap[1:,:]=Xmap
  NZmap[1:,:]=Zmap

  del(Xmap)
  del(Zmap)
  Xmap=NXmap
  Zmap=NZmap

  eq['xmap']=Xmap
  eq['zmap']=Zmap
  eq['jac']=jac
  eq['maxpsi']=sepfrac
  
  if doplot: plot_equilibrium(eq)
  return eq


def plot_equilibrium(eq):

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_aspect('equal')
  fig.set_figheight(4)

  Xmap=eq['xmap'] ; Zmap=eq['zmap']

  #Theta lines
  for i in np.arange(0,len(Xmap[0,:]),4):
    ax.plot(Xmap[:,i],Zmap[:,i])

  #Psi surface
  for i in np.arange(0,len(Xmap[:,0]),5):
    ax.plot(np.append(Xmap[i,:],Xmap[i,0]), np.append(Zmap[i,:],Zmap[i,0]) )
  ax.plot(np.append(Xmap[-1,:],Xmap[-1,0]), np.append(Zmap[-1,:],Zmap[-11,0]) )

  ax.set_title('Surfaces of constant theta and psi (every 5th)');

