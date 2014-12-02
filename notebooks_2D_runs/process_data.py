from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from MITgcmutils import rdmds

def load_UVWT(_Nx,_Nz,_dir,_out_files):
    _U = np.zeros((_Nz,1,_Nx,len(_out_files)));
    _V = np.zeros((_Nz,1,_Nx,len(_out_files)));
    _W = np.zeros((_Nz,1,_Nx,len(_out_files)));
    _T = np.zeros((_Nz,1,_Nx,len(_out_files)));

    for _itr, _out_num in enumerate(_out_files):
        _U[:,:,:,_itr] = rdmds(_dir+'/U',_out_num);
        _V[:,:,:,_itr] = rdmds(_dir+'/V',_out_num);
        _W[:,:,:,_itr] = rdmds(_dir+'/W',_out_num);
        _T[:,:,:,_itr] = rdmds(_dir+'/T',_out_num);

    return np.squeeze(_U),np.squeeze(_V),np.squeeze(_W),np.squeeze(_T);

def plt_field(_depth,_x,_z,_field,_out):
    _Nz = len(_z);
    _total_depth = max(_depth);
    plt.pcolor(_x/1000,_z,np.squeeze(_field[:,:,_out]),cmap='Spectral_r');
    plt.colorbar();
    plt.plot(_x/1000,-_depth-(_total_depth/2)/_Nz,'-k');
    plt.ylim((-_total_depth,0));
