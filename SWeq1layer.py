#################################
# SHALLOW WATER MODEL : 1 LAYER #
#        Dimitri MOREAU         #
#################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from integrator_1_layer import *

#OPTIONAL : PRESENTATION

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['font.size'] = 15

##################################
#GENERAL PARAMETERS

omega = 2*np.pi/(24*3600)
f = 2*omega*np.sin(45)
g = 9.81
cmap = 'RdBu_r'
#f=2e-2 #pour tester coriolis
#H = 1 #replaced by a matrix for topo

integrator = 'LP_m' #'Euler' or 'LP_m'
save_netcf = True #choose if you want to save variables in netCDF

###################################
#PLOT SELECTION

view_1D = False
view_2D = False
view_3D = False

view_cour = False #voir les vecteurs sur view_2D
nb_adim = False #voir les nombres adimensionnels

##################################
#CHOICE OF PERTURBATION

IC = 'eta_detroit' #'decent', 'bidecent', 'cent','eta_detroit'
obstacle = 'none' # 'detroit', 'ile'
rotating_frame = True

###############################
#SPATIAL AND TIME DISCRETISATION
#DECLARATION OF VARIABLES

Lx = 5*1e3
Ly = 5*1e3

Nx = 50
Ny = 50
dx = Lx/Nx
dy = Ly/Ny

#dt = 6 #max euler
dt = 20 #max LP
#time = 1500 #max euler
time = 5000 #max LP

Nt = time/dt
print('Nombre de cellules temporelles : ',Nt) #ne pas dépasser 250 cellules

x = np.arange(0,Lx,dx)
y = np.arange(0,Ly,dy)
t = np.arange(0,time,dt)

u = np.zeros([len(t),len(x),len(y)])
v = np.zeros([len(t),len(x),len(y)])
eta = np.zeros([len(t),len(x),len(y)])

H = np.ones([len(x),len(y)])
#sloped bottom

def gaussian2D(x,y,mx,my):
    from scipy.stats import multivariate_normal
    # Define the parameters of the 2D Gaussian function
   
    sigma_x = np.std(x)//4
    sigma_y = np.std(y)//4
    
    xx,yy=np.meshgrid(x,y)
    pos = np.empty(xx.shape + (2,))
    pos[:, :, 0] = xx
    pos[:, :, 1] = yy
    # Evaluate the 2D Gaussian function at each point in the grid
    rv = multivariate_normal([mx, my], [[sigma_x**2, 0], [0, sigma_y**2]])
    z = rv.pdf(pos)
    print('Max amp : ',np.max(z))

    return z

###################################
#INITIAL CONDITIONS IMPLEMENTATION

corr = 1e6 #correction pour la gaussienne 

if IC == 'decent':
    eta[0,:,:] = gaussian2D(x,y,Lx//6,Ly//6)*corr
elif IC == 'cent':
    eta[0,:,:] = gaussian2D(x,y,Lx//2,Ly//2)*corr
elif IC == 'eta_detroit':
    eta[0,:,:] = gaussian2D(x,y,Lx//2,Ly//6)*corr
elif IC == 'bidecent':
    eta[0,:,:] = gaussian2D(x,y,Lx//2,Ly//2)*corr + gaussian2D(x,y,Lx//6,Ly//6)*corr

########################################################
#COMPUTATION

if integrator == 'Euler':
    if rotating_frame == False:
        f = 0
        print('Coriolis parameter : ',f)
        u,v,eta = get_SW_euler_R(t,x,y,Ny,u,v,eta,dx,dy,dt,g,H,f,obstacle)
    elif rotating_frame == True:
        print('Coriolis parameter : ',f)
        u,v,eta = get_SW_euler_R(t,x,y,Ny,u,v,eta,dx,dy,dt,g,H,f,obstacle)

elif integrator == 'LP_m':
    if rotating_frame == False:
        f = 0
        print('Coriolis parameter : ',f)
        u,v,eta = get_SW_LP_mixte_R(t,x,y,Ny,Nt,u,v,eta,dx,dy,dt,g,H,f,obstacle)
    elif rotating_frame == True:
        print('Coriolis parameter : ',f)
        u,v,eta = get_SW_LP_mixte_R(t,x,y,Ny,Nt,u,v,eta,dx,dy,dt,g,H,f,obstacle)
  
###############################################################
#PLOT

if view_1D == True:
    plt.figure()
    n = 25
    eps=1
    for i in range(len(t)):
        plt.plot(x,eta[i,:,n]+H[:,n])
        plt.xlabel(r'$x$')
        plt.ylabel(r'$H$')
        plt.ylim(np.min(eta),np.max(eta)+eps)
        plt.title(str(i)+'/'+str(len(t)-1)+r' at $y=$'+str(n))
        plt.pause(0.01)
        plt.clf()
    plt.plot(x,eta[-1,:,n]+H[:,n])
    plt.xlabel(r'$x$')
    plt.ylabel(r'$H$')
    plt.title(str(i)+'/'+str(len(t)-1)+r' at $x=$'+str(n))


if view_2D == True:
    plt.figure()
    for i in range(len(t)):
        plt.pcolormesh(x,y,eta[i,:,:],cmap=cmap,vmin=np.min(eta),vmax=np.max(eta))
        #plt.contourf(x,y,eta[i,:,:],cmap=cmap,vmin=np.min(eta),vmax=np.max(eta))
        plt.colorbar(extend="both",label=r'$\eta$')
        plt.xlabel(r'$x$')
        plt.ylabel(r'$y$')
        plt.title(str(i)+'/'+str(len(t)-1))
        #plt.savefig('/home/dmoreau/Bureau/PDE/gif_shallow/im_'+str(i+1))
        plt.pause(0.01)
        plt.clf()
    plt.pcolormesh(x,y,eta[-1,:,:],cmap=cmap,vmin=np.min(eta),vmax=np.max(eta))


if view_3D == True:
    fig, ax = plt.subplots(1,1,figsize=(10,7),subplot_kw={"projection": "3d"})
    xx,yy=np.meshgrid(x,y)
    for i in range(len(t)):

        fig.suptitle(str(i)+'/'+str(len(t)-1))
        fig3=ax.plot_surface(xx,yy,eta[i,:,:]+H,cmap=cmap)
        ax.set_zlim(np.min(eta),np.max(eta))

        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        ax.set_zlabel('$H$')


        plt.pause(0.01)
        ax.clear()
    fig3=ax.plot_surface(xx,yy,eta[-1,:,:]+H,cmap=cmap)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_zlabel('$H$')
    ax.set_zlim(np.min(eta),np.max(eta))

if nb_adim == True:
    U=u
    Ro = U/(f*Lx)
    Fr = U/(np.sqrt(g*H))
    print(Ro)
    fig,(ax) = plt.subplots(1,2,figsize=(15,7))
    for i in range(len(t)):
        fig.suptitle(str(i)+'/'+str(len(t)-1))

        fig1=ax[0].pcolormesh(x,y,Ro[i,:,:],cmap=cmap,vmin=np.min(Ro),vmax=np.max(Ro))
        fig2=ax[1].pcolormesh(x,y,Fr[i,:,:],cmap=cmap,vmin=np.min(Fr),vmax=np.max(Fr))

        ax[0].set_title(r'$R_o = \frac{u}{f.L}$')
        ax[1].set_title(r'$F_r = \frac{u}{\sqrt{g.H}}$')

        plt.pause(0.01)
        ax[0].clear()
        ax[1].clear()


if view_cour == True:
    plt.figure()
    for i in range(len(t)):
        plt.pcolormesh(x,y,eta[i,:,:],cmap=cmap,vmin=np.min(eta),vmax=np.max(eta))
        plt.colorbar(extend="both",label=r'$\eta$')
        plt.quiver(x,y,u[i,:,:],v[i,:,:])
        plt.xlabel(r'$x$')
        plt.ylabel(r'$y$')
        plt.title(str(i)+'/'+str(len(t)))

        plt.pause(0.01)
        plt.clf()
    
plt.show()

def save_nc(u,v,eta):
    from netCDF4 import Dataset
    # Specify the netCDF file name
    netcdf_filename = 'output.nc'
    # Create a netCDF file in write mode
    with Dataset(netcdf_filename, 'w', format='NETCDF4') as ncfile:
        # Create dimensions matching the shape of your arrays ## common for all variables in t,x,y
        dim11 = ncfile.createDimension('time', eta.shape[0])
        dim12 = ncfile.createDimension('x', eta.shape[1])
        dim13 = ncfile.createDimension('y', eta.shape[2])
        # Create variables and store the arrays
        var1 = ncfile.createVariable('eta', 'f8', ('time','x', 'y'))
        var1[:] = eta
        var2 = ncfile.createVariable('u', 'f8', ('time','x', 'y'))
        var2[:] = u   
        var3 = ncfile.createVariable('v', 'f8', ('time','x', 'y'))
        var3[:] = v      

if save_netcf == True:
    save_nc(u,v,eta)


'''
#################################
#BC : domaine fermé

#rigide : 
#u nul en 0 et -1 x
#v nul en 0 et -1 y

#rigide et visqueux :
#u nul en 0 et -1 x et 0 et -1 y
#v nul en 0 et -1 y et 0 et -1 x

########################################################
#COMPUTATION'''