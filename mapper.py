#mapper program. Expects eq variable to be defined with eqdsk function.

import scipy.integrate as integrate
import scipy.fftpack as sft
from scipy import integrate
import scipy.interpolate

mu0=4.*np.pi*1.e-7
R=eq.get('r')
Z=eq.get('z')
B,gpsi,fRZ,Rv,Zv=eqdsk.getModB(eq)
psi=eq.get('psizr').T
curtor=[]
area=[]

ifrhopol=True #use root psipol mesh instead of psipol (eg for torlh)

def find_cut(x,y, rmaxis, zmaxis):
  #from S. Shiraiwai to find crossing going counter clockwise
  for k in range(len(y)-1):
       km = k-1
       if y[km] < zmaxis and y[k] > zmaxis:
         return k
    #('findcut',x,y)
  return -1

#ensure periodicity with fft and reduce noise, esp near axis

# get flux surface on fine rectangular mesh
newR=np.arange(min(R),max(R),0.01)
newZ=np.arange(min(Z),max(Z),0.01)

RR,ZZ=np.mgrid [min(R):max(R):129j, min(Z):max(Z):129j ] #129x129 2D coordinate meshes
RR_finer,ZZ_finer=np.mgrid [min(R):max(R):200j, min(Z):max(Z):200j ] #200x200 2D coordinate meshes

spline_psi = scipy.interpolate.RectBivariateSpline(R,Z,psi.T,bbox=[np.min(R),np.max(R),np.min(Z),np.max(Z)],kx=5,ky=5)
psi_int=spline_psi.ev(RR_finer,ZZ_finer)
psi_int_r=spline_psi.ev(RR_finer,ZZ_finer,dx=1)
psi_int_z=spline_psi.ev(RR_finer,ZZ_finer,dy=1)
grad_psi=np.sqrt(psi_int_z**2+psi_int_r**2)

#check for save
import pickle
try:
  f=open(rundir+eqfileobj, "rb")
except FileNotFoundError:
  print('Recalculating mapping')

  #generated mapped mesh size:
  npsi=128
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

  nidx=100
#reference psi
  fity=np.linspace(eq['simag'],eq['sibry'],nidx)

  if ifrhopol:
    rhopol = np.linspace(np.sqrt(np.abs(eq['simag'])),np.sqrt(np.abs(eq['sibry'])),nidx)*np.sign(fity)
    fity = rhopol**2*np.sign(fity) #reference psipol consistent with uniform rhopol
    #rhopol is linear so just declare it
    rhopol = np.linspace(0,1,npsi) 
    #rhopol = (rhopol-np.min(rhopol))/(np.max(rhopol)-np.min(rhopol)) #normalize
    eq['rhopolmap']=rhopol  #sqrt norm rho pol for map size npsi, linear spaced

  fitx = np.linspace(0,nidx,nidx)
#get psimesh
  fitxvals = np.linspace(0, nidx, npsi)
  psimesh = np.interp(fitxvals, fitx, fity) #this just stretches fity to size of npsi

  rmaxis = eq['rmaxis']
  zmaxis = eq['zmaxis']
  eq['psipolmap'] = psimesh

  c_pprime = np.interp( psimesh, eq['fluxGrid'], eq['pprime'] )
  c_ffprime = np.interp( psimesh, eq['fluxGrid'], eq['ffprim'] )

#Extract contours and values for flux coordinate system.
#contours go counter-clockwise, which we want
#contours don't necessarily start at y=0., so rebase
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_aspect('equal')
  fig.set_figheight(6)
  psi_cs=plt.contour(RR_finer,ZZ_finer,psi_int,levels=psimesh)

#Get Psi contours, careful to exclude field coils
  psixy=[]
  for p in psi_cs.collections:
      for pp in p.get_paths():
        v = pp.vertices
        x = v[:,0]
        y = v[:,1]
        if np.abs(np.average(y))<0.20*np.max(eq['z']): #only keep core plasma contours
          psixy.append( (x,y) )

#Define uniform theta mesh
  uni_theta=np.linspace(0,2.0*np.pi,ntheta,endpoint=False)

#Set up X(psi,theta) Y(psi,theta) and initialize with magnetic axis point
#actually, do NOT add axis point
  eq_x=[]
#eq_x.append(np.array([rmaxis]*ntheta))
  eq_y=[]
#eq_y.append(np.array([zmaxis]*ntheta))

  points = np.array( (RR_finer.flatten(), ZZ_finer.flatten()) ).T
  psivalues = grad_psi.flatten()
  Bpoints = np.array( (RR.flatten(), ZZ.flatten()) ).T
  Bvalues   = B.flatten()

  from scipy.interpolate import griddata
  print('points',points.shape,psivalues.shape)
  c_idx=-1
  for cx,cy in psixy:
      c_idx+=1
      print('cx',c_idx, len(cx),len(psixy)) #just to see progress
      #for each surface, low pass filter to central 8+ DC Fourier modes
      #remove last element for fft since it is equal to the first element
      #size of (cx,cy) is  variable
      minmodes = 8
      maxmodes = 20
      nmodes=max( minmodes,int(len(cx)/np.float(maxmodes) ))
      fftx=sft.fft(cx[0:-1])
      fftx[int(nmodes/2)+1:-int(nmodes/2)]=0
      filtered_cx=sft.ifft(fftx).real
        
      ffty=sft.fft(cy[0:-1]) #fft makes filtered_cx/y complex
      ffty[int(nmodes/2)+1:-int(nmodes/2)]=0
      filtered_cy=sft.ifft(ffty).real

    #shift contours and insert outer midplane value so we have theta=0 value
    #an important step to remove noise

      idx=find_cut(filtered_cx,filtered_cy,rmaxis,zmaxis)
      if (idx>=0):
          filtered_cx=np.roll(filtered_cx,-idx)
          filtered_cy=np.roll(filtered_cy,-idx)
          #print('rolled',filtered_cy[0],filtered_cy[-1])
    #interpolate |grad psi| and B onto this surface
    #this significantly slows down this routine. I think this can be moved outside of the 
    #loop if cx,cy are saved.
      c_gradpsi = griddata( points, psivalues, (filtered_cx,filtered_cy), method='cubic' )
      c_B       = griddata( Bpoints,   Bvalues, (filtered_cx,filtered_cy), method='cubic' )

    #derivative from fft needs factor of 2pi
      df_dx=sft.diff(filtered_cx)*np.pi*2.0/len(filtered_cx)
      df_dy=sft.diff(filtered_cy)*np.pi*2.0/len(filtered_cy)
      dl=np.sqrt(df_dx**2+df_dy**2) #These two steps could be done with FFT too.

      #flux surface integrals here
      print(c_idx, np.sum(dl) , integrate.simps(dl), np.average(1/c_gradpsi)*(psimesh[c_idx]-psimesh[c_idx+1]))
      c_curtor = -integrate.simps(
          dl*filtered_cx/c_gradpsi * ( c_pprime[c_idx] + c_ffprime[c_idx]/filtered_cx**2/mu0 ) 
      )
      c_area = integrate.simps(
          dl/c_gradpsi
      )


      #if straight field line, multiply dl by 1/R*|gradpsi|
      dtheta=dl/np.abs(c_gradpsi*filtered_cx)
      L=integrate.cumtrapz(dtheta,initial=0)/len(dtheta) 
      #L=np.insert(L,0,0)

      #now put each on same theta mesh, Jacobian selection
      this_theta=L/np.max(L)*2.0*np.pi
      t_map=scipy.interpolate.interp1d(this_theta,filtered_cx,kind='cubic')
      eq_x.append(t_map(uni_theta))
      t_map=scipy.interpolate.interp1d(this_theta,filtered_cy,kind='cubic')
      eq_y.append(t_map(uni_theta))
      curtor.append (c_curtor)
      area.append( c_area )


  #add origin term
  area.insert(0,0.)          #area of origin is 0. 
  curtor.insert(0,curtor[0]) #current density maximum at origin
  curtor=np.array(curtor)
  area=np.array(area)
  Xmap=np.real(np.vstack(eq_x))
  Zmap=np.real(np.vstack(eq_y))
  eq['darea_dpsi']=area
  eq['dI_dpsi']=curtor
  #add origin
  NXmap=np.zeros([npsi,ntheta])
  NZmap=np.zeros([npsi,ntheta])
  NXmap[0,:]=eq['rmaxis']
  NZmap[0,:]=eq['zmaxis']
  NXmap[1:,:]=Xmap
  NZmap[1:,:]=Zmap
  del(Xmap)
  del(Zmap)
  Xmap=NXmap
  Zmap=NZmap

  print('shapes of mapped arrays: ',Xmap.shape,Zmap.shape)
  #save this in pickle file
  eq['xmap']=Xmap
  eq['zmap']=Zmap
  

  import pickle  
  with open(rundir+eqfileobj, "wb")  as file_pi:
    pickle.dump([Xmap,Zmap,psixy,curtor,eq], file_pi);

else:
    try:
      Xmap,Zmap,psixy,curtor,eq = pickle.load(f)
      print('loaded previous mapping')
    finally:
        f.close()
