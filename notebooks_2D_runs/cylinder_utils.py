import numpy as np
import scipy.integrate as integrate

def find_dvtheta_dr(x,z,u,w,ux,uz,wx,wz):
    ''' find_dvtheta_dr
    Takes:
      x: x position value (float)
      z: z position value
      u: u velocity value (float)
      w: w velocity value
      u_x: x derivative of u value
      u_z: z derivative of u value
      w_x: x derivative of w value
      w_z: z derivative of w value
    Returns:
      dvtheta_dr: floating point value of the polar derivative

    This function takes the Cartesian coordinates and derivatives of the u and w velocity fields and
    converts them into the polar coordinate derivative described physically by:
    the change in along cylinder velocity (vtheta) with respect to distance from the cylinder (r)
    '''
    return -(-wx*x**2+x*w-x*z*wz+x*z*ux-z*u+z**2*uz)/(x**2+z**2)**(1.5)

def find_bdy(X,Z,R0,T):
    ''' find_bdy
    Takes:
      X: 1D array of x position values
      Z: 1D array of z position values
      R0: physical size of cylinder (2.5cm)
      T: the temparture field; used to determine boundary
    Returns:
      bd_pts: 1D array of size 2 array values specifying the [z index, x index] values of boundary
              points.

    This function scans the temperature field for points that are 0 degrees C. The
    temperature field is not stepped forward in time, and so points that are 10 degrees C
    are interior points and points that are 0 degrees C are boundary points. A single boundary
    point is specified as the x index and z index of an interior cell that is adjacent to a
    boundary cell.    
    '''
    bd_pts = []

    Nx = X.shape[0]
    Nz = Z.shape[0]    

    for itr,ii in enumerate(X):
        zind = Nz-1
        for itr2,jj in enumerate(Z):
            if T[itr2,itr] == 0:
                zind = itr2-1
                break
        if X[itr] > -R0 and X[itr] < R0:
            bd_pts.append([zind,itr])
            
    return bd_pts

def calc_dvtheta_dr(X,Z,U,W,R0,Nx,Nz,bd_pts):
    ''' calc_dvtheta_dr
    Takes:
      X: 1D array of x position values
      Z: 1D array of z position values
      U: 2D array of u values (ordered [Z,X])
      W: 2D array of w values
      R0: physical radius of cylinder (2.5cm)
      Nx: number of x points in domain
      Nz: number of z points in domain
      bd_pts: the boundary points, as calculated in find_bdy
    Returns:
      theta_bd: the theta value of boundary points (in degrees, measured from the left
                stagnation point)
      dvtheta_dr: the derivative at the boundary points
 
    This function calculates the value of dvtheta_dr for each boundary cell.
    '''
    dx = (X[-1]-X[0])/Nx
    dz = abs((Z[-1]-Z[0])/Nz)

    dvtheta_dr_bd = []
    theta_bd = []

    for ii in bd_pts:
        if X[ii[1]] > -R0 and X[ii[1]] < R0:
            xind = ii[1]
            zind = ii[0]
            xx = X[ii[1]]
            zz = Z[ii[0]]
            u_bd = U[zind,xind]
            w_bd = W[zind,xind]
            if ii < 0:
                ux_bd = (U[zind,xind]-U[zind,xind-1])/dx
                wx_bd = (W[zind,xind]-W[zind,xind-1])/dx
            else:
                ux_bd = (U[zind,xind+1]-U[zind,xind])/dx
                wx_bd = (W[zind,xind+1]-W[zind,xind])/dx
            uz_bd = (U[zind,xind])/dz
            wz_bd = (W[zind,xind])/dz
            dvtheta_dr_bd.append(find_dvtheta_dr(xx,zz,0,0,ux_bd,uz_bd,wx_bd,wz_bd))
            theta_bd.append(np.arccos(-xx/np.sqrt(xx**2+zz**2))*360.0/(2*np.pi))
            
    return theta_bd,dvtheta_dr_bd

def find_zero(theta_bd,dvtheta_dr_bd):
    ''' find_zero
    Takes:
      theta_bd: 1D array of theta values at boundary
      dvtheta_dr_bd: 1D array of derivative at boundary
    Returns:
      theta_0: the value of theta where the derivative is zero

    This function approximately determines where the zero of the derivative at the boundary occurs.
    Effectively this finds the separation point.
    '''
    prev_val = dvtheta_dr_bd[0]
    curr_val = -1
    theta_0 = 180
    for ii in range(1,len(dvtheta_dr_bd)):
        curr_val = dvtheta_dr_bd[ii]
        if prev_val < 0 and curr_val > 0:
            print "Found zero!"
            theta_0 = theta_bd[ii-1]
            break
        prev_val = curr_val
    return theta_0

def calc_phi(X,Z,U,W,bd_pts):
    ''' calc_phi
    Takes:
      X: 1D array of X values
      Z: 1D array of Z values
      U: 2D array of U values in domain (ordering: [Z,X])
      W: 2D array of W values in domain (ordering: [Z,X])
      bd_pts: UNUSED, specifies boundary points
    Returns:
      phi: 2D array of the velocity potential, phi
    
    This function numerically integrates the U and W fields to form the phi field.
    Method of integration is the cumulative trapezoidal method implemented in scipy.integrate
    package. The u and w fields are extended based on theoretical data (U = 0.01, W = 0) so
    that the integration returns fields of the same size. The mean is subtracted off of the
    integrated fields (because the potential function is unique up to an undetermined constant).
    '''
    
    Nz,Nx = U.shape;
    
    u_extx = np.append(U,np.zeros((Nz,1))+0.01,axis=1)
    w_extz = np.append(np.zeros((1,Nx)),W,axis=0)
    
    phix = integrate.cumtrapz(u_extx,np.append(X,X[-1]+(X[1]-X[0])),axis=1)
    phiz = integrate.cumtrapz(w_extz,np.append(Z,Z[-1]+(Z[1]-Z[0])),axis=0)
    
    phi = phix + phiz
    
    phi = phi - np.mean(phi)
    
    return phi
