import numpy as np

def find_dvtheta_dr(x,z,u,w,ux,uz,wx,wz):
    return -(-wx*x**2+x*w-x*z*wz+x*z*ux-z*u+z**2*uz)/(x**2+z**2)**(1.5)

def find_bdy(X,Z,R0,T):
    bd_pts = []
    
    for itr,ii in enumerate(X):
        zind = 127
        for itr2,jj in enumerate(Z):
            if T[itr2,itr] == 0:
                zind = itr2-1
                break
        if X[itr] > -R0 and X[itr] < R0:
            bd_pts.append([zind,itr])
            
    return bd_pts

def calc_dvtheta_dr(X,Z,U,W,R0,Nx,Nz,bd_pts):
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
    prev_val = dvtheta_dr_bd[0]
    curr_val = -1
    theta_0 = 180
    for ii in range(1,len(dvtheta_dr_bd)):
        curr_val = dvtheta_dr_bd[ii]
        if prev_val < 0 and curr_val > 0:
            print "Found zero!"
            theta_0 = (theta_bd[ii]+theta_bd[ii-1])/2
            break
        prev_val = curr_val
    return theta_0
