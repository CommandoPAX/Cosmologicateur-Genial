#################Attention la BoxSize soit etr emodifie a deux endroits (dont un hard code un peu cache)###################
import yt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import interpolate
import MAS_library as MASL
import h5py
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy.stats import binned_statistic

from mpl_toolkits.mplot3d import Axes3D
import warnings

plt.rcParams['figure.figsize'] = [20, 7]


matplotlib.rcParams.update({'font.size': 20})

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams["figure.facecolor"]='w'
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

BoxSize  = 500
grid     = 32                   #grid size
ptypes   = [1]                   #CDM + neutrinos
MAS      = 'CIC'                   #Cloud-in-Cell
do_RSD   = False                   #dont do redshif-space distortions
axis     = 0                       #axis along which place RSD; not used here
verbose  = False   #whether print information on the progress

#Dp = 0.025223
N=32


class Params():
    '''Class of parameters'''
    def __init__(self):
    
        self.h=0.67742
        self.omega_b=0.048891054
        self.omega_cdm=0.261008
        self.omega_ncdm= 0 #0.00140718
        self.N_ncdm=0
        self.m_ncdm=0.0
        
        self.omega_r= 0 #9.16714e-05
        self.omega_k=0
        self.z=50

        self.A_s = 2.1064e-9
        self.n_s = 0.96822
        self.k_pivot = 0.05 # 1/Mpc
        
        self.N=32                  # grid size                                
        self.kmin=2*np.pi/BoxSize

params=Params()

k_min= 2*np.pi/BoxSize
k=np.linspace(-(N//2)*k_min,N//2*k_min,N+1,dtype=np.float64)
k_grid=np.array(np.meshgrid(k,k,k,sparse=True,indexing='ij'),dtype=object) 
#kmod=np.sqrt(k_grid[0][params.N//2:]**2+k_grid[1]**2+k_grid[2]**2)


def ifft(field):
    field[0,params.N//2,params.N//2]=0
    return np.fft.irfftn(np.fft.ifftshift(field.transpose()[:-1,:-1],axes=(0,1)),(params.N,params.N,params.N) )

def fft(f_field): # fast fourier transform
    field=np.zeros((params.N//2+1,params.N+1,params.N+1),dtype=complex)
    field[:,:-1,:-1]=np.fft.fftshift(np.fft.rfftn(f_field),axes=(0,1)).transpose()
    field[:,-1],field[:,:,-1]=field[:,0],field[:,:,0]
    return field

def W(grid):
    l = BoxSize/params.N
    k1,k2,k3=grid[0][params.N//2:]*l/2/np.pi,grid[1]*l/2/np.pi,grid[2]*l/2/np.pi
    return (np.sinc(k1)*np.sinc(k2)*np.sinc(k3))**2

def WW(grid):
    R=8
    Rkmod=np.sqrt(grid[0][params.N//2:]**2+grid[1]**2+grid[2]**2)*R
    return 3 * (np.sin(Rkmod) - Rkmod*np.cos(Rkmod))/Rkmod**3



def sigma8(field,k_grid):
    fdelta = fft(field)

    Wcic = W(k_grid)
    Wth = WW(k_grid)

    fdelta = fdelta * Wth / Wcic

    delta = ifft(fdelta)
    return np.sqrt(np.mean(delta**2))
#print( r'$\sigma_8 = $'+str(np.sqrt(np.mean(delta**2)) ))

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    


<<<<<<< HEAD
    index = 3
    DATA = yt.load("../2024-03-12 15:55:46 - WMDPGN32/output_00003/info_00003.txt")
=======
    index = 2
    DATA = yt.load("./RESULT/2024-03-12 15:55:46 - WMDPGN32/output_00003/info_00003.txt")
>>>>>>> aae99dc0d20c674542335154c92f0fa0151fe580
    DATA.all_data().to_dataframe(["particle_position_x","particle_position_y","particle_position_z","particle_mass"])

    grid = 32    #grid size
    pBoxSize = DATA.domain_width.in_units('Mpccm/h') #Mpc/h
    BoxSize = pBoxSize[0].value #Mpc/h
    Rayleigh_sampling = 1     #whether sampling the Rayleigh distribution for modes amplitudes
    threads = 1      #number of openmp threads
    verbose = False   #whether to print some information
    axis = 0
    MAS = 'CIC'

    ad=DATA.all_data()
    pos = ad['particle_position'].astype(np.float32)*BoxSize

    # define 3D density fields
    delta = np.zeros((grid,grid,grid), dtype=np.float32)

    # construct 3D density field
    MASL.MA(pos.astype(np.float32), delta, BoxSize, MAS, verbose=verbose)

    # at this point, delta contains the effective number of particles in each voxel
    # now compute overdensity and density constrast
    delta /= np.mean(delta, dtype=np.float64);  delta -= 1.0
            
        
    res=sigma8(delta,k_grid)
#np.savetxt("S8.dat",res)
print("sigma8=", res)









