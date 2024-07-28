#!/usr/bin/env python

"""
Python class for reading in cql3d output files. Returns handle to cql3d instance.
"""
#for math and numerical functions
import pylab as m
import matplotlib.pyplot as plt
#netcdf interface
import scipy.io.netcdf as nc
import numpy as np

#UPdated May2011 for scipy 0.10 with .data netcdf attribute
class cql3d:
    "I/O and plotting routines for the CQL3D code"

    cqlname='undefined'
    cqlrfname='undefined'
    cqlhdl=0
    cqlrfhdl=0
    cqldict={}

    def __init__(self, cqlroot):
        self.cqlname=cqlroot+'.nc'
        self.cqlrfname=cqlroot+'_rf.nc'
##Open the cql3d netcdf file read only
        try:
            self.cqlhdl = nc.netcdf_file(self.cqlname,'r')
        except IOError:
            print ('CRITICAL',self.cqlname,'not found.')
            self.cqlhdl = -1

        try:
            self.cqlrfhdl = nc.netcdf_file(self.cqlrfname,'r')
        except IOError:
            print (self.cqlrfname,'not found, relevant fns disabled.')
            self.cqlrfhdl = -1
        return
    
    def info( self ):
        "Prints a list of the contents of present CQL3D output files"

        if (self.cqlhdl != -1):
            print ('The cql file, ',self.cqlname,', contains:')
            print ('----------------------------------')
            print ("The global attributes: ",self.cqlhdl.dimensions.keys())
            print ("File contains the variables: ", self.cqlhdl.variables.keys())

        if (self.cqlrfhdl != -1):
            print ('The cql file, ',self.cqlrfname,', contains:')
            print ('----------------------------------')
            print ("The global attributes: ",self.cqlrfhdl.dimensions.keys())
            print ("File contains the variables: ", self.cqlrfhdl.variables.keys())

        return


    def cqlparam( self ):
        "Fill a dictionary of cql scalar values, not implemented."
        self.cqldict['comment']='Dictionary of cql3d outputs'
        self.cqldict['rmaj']=cqltest.cqlhdl.variables['radmaj'].getValue()
        self.cqldict['rmag']=cqltest.cqlhdl.variables['rmag'].getValue()

        return



#option of using radius index of value (nearest). Time step needed?
    def fplot( self, irad=-2, itime=-1, fig='none',var='none',ptype='contour' ):
        """Contour plot of the distribution function at given radius.
        irad may be either the radial index or value. Also plots other
        3D velocity space variables. Use type='line' to get 8 pitch angle
        slices""" 

        figscale=3 #scale up figure (1 in height originally)

        if var=='none':
            f     = self.cqlhdl.variables['f']

        if var=='B':
            f=self.cqlrfhdl.variables['rdcb']

        if var=='rayB':
            f=self.cqlrfhdl.variables['urfb']


        u     = self.cqlhdl.variables['x'][:]
        pitch = self.cqlhdl.variables['y'][:]
        rya   = self.cqlhdl.variables['rya'][:]

#f is on a polar coordinate mesh. We need to generate a 2D cartesian mapping
#to vpar=u cos(pitch[ir,:]) and vperp = u sin(pitch[ir,:])
#the numpy meshgrid command does this for us.

#just use the last surface, we are assuming theta mesh is the same on each surface
        nr = self.cqlhdl.dimensions['rdim']
        if (irad == -2):
            irad = nr/2

        if (type(irad)==float):
            irad=int(np.interp(irad,rya[:],np.arange(np.size(rya))))

        if (irad > nr-1):
            print ("irad too large, using nr/2=0.5",irad)
            irad = nr/2

        r,t = np.meshgrid(u,pitch[irad,:])
        vpar  = np.transpose(r*np.cos(t))
        vperp = np.transpose(r*np.sin(t))

        if (ptype=='line'):
            fig=plt.figure()
            nt = self.cqlhdl.dimensions['ydim']
            c  = 2.99792458e10
            vnorm = self.cqlhdl.variables['vnorm'].getValue()
            uscaled = np.array(u)*vnorm/c
            for iy in np.arange(0,nt,nt/8):
                plt.plot(uscaled,np.log10(f[irad,:,iy]),label=str(int(iy/float(nt)*8))+r'/8 $\pi$')
                
            plt.axis( ymin=0 )
            plt.title(str(f.long_name)+" ["+str(f.units)+"]",size=24)
            plt.xlabel(r'$\gamma v_{||}/c$',size=22)
            plt.ylabel(r'$\gamma v_\bot/c$',size=22)
            plt.legend(loc=1)
            plt.draw()
            
            return


#for half plane plot aspect 2x1 and add space for colorbar legend
        if (fig=='none'):
            fig=plt.figure(figsize=(2.8*figscale,1*figscale)) 
            fig.subplots_adjust(bottom=0.15,hspace=0.3)
#else type(fig)==matplotlib.figure.Figure then ok
#could also use gcf() to get current figure.

        ax = fig.add_subplot(111) #plot on existing canvas
        v=(np.arange(30)/3.+7) #contouring levels, 6 appropriate for /cc
        csf = ax.contourf(vpar,vperp,np.log10(abs(f[irad,:,:])+1))
#        cs = ax.contour(vpar,vperp,f[irad,:,:],v)
# make sure aspect ratio preserved
        ax.set_aspect('equal')

# draw a circle around the edge of the plot.
        rmax = max(u)
        theta=np.arange(100)*np.pi/99.
        ax.plot(rmax*np.cos(theta),rmax*np.sin(theta),'k')
        plt.title(str(f.long_name)+" ["+str(f.units)+"] at r="+str(rya[irad])[0:4],size=24)
        plt.xlabel(r'$\gamma v_{||}/c$',size=22)
        plt.ylabel(r'$\gamma v_\bot/c$',size=22)
        plt.axis([-1,1,0,1])
#        cax = plt.axes([0.85, 0.2, 0.02, 0.6])
        cbar=plt.colorbar(csf,format='%4.2f',ax=ax) #,fraction=0.05
        cbar.ax.set_ylabel('log10 levels')
        plt.draw()
#,orientation='horizontal')

        return
    

    def pltpower( self, ispecies=0, itime=-1 ):
        "Plot power deposition."

        rfpwr = self.cqlhdl.variables['rfpwr']
        ptitle = rfpwr.long_name0
        plabel = rfpwr.units
        if (ispecies == 4):
            ptitle = rfpwr.long_name1
            plabel = 'Watts' #not in attributes
        self.__rhoplot( rfpwr[itime,ispecies,:], ptitle, plabel )

        return

    def __rhoplot( self, pvar, ptitle='', plabel='' ):
        "Hidden internal 1D plot vs rho"

        rya=self.cqlhdl.variables['rya']
        plt.plot(rya, pvar, 'g-o')
        plt.title(ptitle)
        plt.ylabel(plabel)
        plt.xlabel(rya.long_name)

        return


    def __getir( self, irad ):
        "Gets radial index and value, given one of them, not implemented."

        return

#interfaces to ray data start here.
    def pltrays( self, pctpwr=0.001, iskip=10, fig=False):
        from matplotlib.collections import LineCollection
        """Plots rays. 
        pctpwr: gives power percent cutoff, defaults to nearly 100%
        iskip: plots every iskip'th ray
        """

        wr=self.cqlrfhdl.variables['wr']
        wz=self.cqlrfhdl.variables['wz']
        delpwr=self.cqlrfhdl.variables['delpwr']
#number of rays, also the size of the first index in wr,wz,and delpwr
        nrays=self.cqlrfhdl.dimensions['nrays']
#number of ray steps, also the size of the second index in wr,wz,and delpwr
        neltmax=self.cqlrfhdl.dimensions['neltmax']
#e.g. wz.dimensions = ('nrays', 'neltmax')
        rmag=self.cqlhdl.variables['rmag'].getValue()

        if not fig:
            fig, ax = plt.subplots()
        else:
            ax=plt.gca()
        
        for iray in range(0,nrays,iskip):
            lastidx = min(np.concatenate( (np.nonzero(delpwr[iray,:] < 
                          pctpwr*delpwr[iray,0])[0],[neltmax-1]) ))
            cwz=0.0
            cwr=rmag
            points=np.array([wr[iray,:lastidx]-cwr,wz[iray,:lastidx]-cwz]).T.reshape(-1,1,2)
            segments=np.concatenate( [points[:-1], points[1:]], axis=1 )
            
            lc = LineCollection(segments, cmap=plt.cm.rainbow)
#                                norm=plt.Normalize(0, pmax))
            lc.set_array(delpwr[iray,:lastidx])
            lc.set_linewidth(1)
            ax.add_collection(lc)
#            plt.plot(wr[iray,:lastidx]-rmag,wz[iray,:lastidx])
        fig.colorbar(lc)
        ax.autoscale()
        plt.axes().set_aspect('equal')
        plt.show()

        return


    
if __name__ == '__main__':
    from pylab import *
    from matplotlib import *
    from  matplotlib.pyplot import *
    import cql3d_mod

    cqltest = cql3d_mod.cql3d('cqlout1060728transp.0.0')

##Read some scalar values using netcdf interface
    rmaj=cqltest.cqlhdl.variables['radmaj'].getValue()
    rmag=cqltest.cqlhdl.variables['rmag'].getValue()

    print ("The major and magnetic axis are:",rmaj,rmag)

##Easier way using cql3d dictionary:


##Read some arrays and convert them to python numpy representation
    rfpower=cqltest.cqlhdl.variables['rfpwr']
    print ("The dimension names of rfpower are:",rfpower.dimensions)
    print ("and its shape is:",rfpower.shape)

#Example of iterating over dimensions of variable and getting more info:
    for dimname in rfpower.dimensions:
        print ("The size of dimension",dimname,"is :",cqltest.cqlhdl.dimensions[dimname])

#file info
    cqltest.info()

#turn interactivity on
    ion()
#power plot, first power density, then integrated power
    figure()
    cqltest.pltpower(0)
    figure()
    cqltest.pltpower(4)

#contour of distribution function, irad can be the radial index or normalized radius
    f=figure()
    cqltest.fcontour(irad=0.5,fig=f)

#plot some rays
    figure()
    cqltest.pltrays()

    show()



