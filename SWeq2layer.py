##################################
# SHALLOW WATER MODEL : 2 LAYERS #
#         Dimitri MOREAU         #
##################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from integrator_2_layers import *

#OPTIONAL : PRESENTATION

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['font.size'] = 15

##################################
#GENERAL PARAMETERS

omega = 2*np.pi/(24*3600)
f = 2*omega*np.sin(45)
g = 9.81
cmap='RdBu_r'
rho2 = 1020
rho1 = 1010
rho0 = 1013
g2= g*(rho1 - rho2)/rho0 #reduced gravity
print('Reduced gravity : ',g2)
#H1=1 #replaced with a matrix for topo
#H2=1
#f=1e-2 #pour tester

integrator = 'LP' #'Euler' or 'LP' or 'LP_m'
save_netcf = True #choose if you want to save variables in netCDF

###################################
#PLOT SELECTION

view_1D = False
view_2D = False
view_3D = True

view_cour = False

##################################
#CHOICE OF PERTURBATION

IC = 'decent' #'decent', 'bidecent', 'cent','eta_detroit'
obstacle = 'none' #'detroit', 'ile'
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
print('Nombre de cellules temporelles : ',Nt)

x = np.arange(0,Lx,dx)
y = np.arange(0,Ly,dy)
t = np.arange(0,time,dt)

u1 = np.zeros([len(t),len(x),len(y)])
v1 = np.zeros([len(t),len(x),len(y)])
eta1 = np.zeros([len(t),len(x),len(y)])

u2 = np.zeros([len(t),len(x),len(y)])
v2 = np.zeros([len(t),len(x),len(y)])
eta2 = np.zeros([len(t),len(x),len(y)])

detadt = np.zeros([len(t),len(x),len(y)])

H1 = np.ones([len(x),len(y)])*1
H2 = np.ones([len(x),len(y)])*1

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
    eta1[0,:,:] = -gaussian2D(x,y,Lx//6,Ly//6)*corr
elif IC == 'cent':
    eta1[0,:,:] = -gaussian2D(x,y,Lx//2,Ly//2)*corr
elif IC == 'eta_detroit':
    eta1[0,:,:] = -gaussian2D(x,y,Lx//2,Ly//6)*corr
elif IC == 'bidecent':
    eta1[0,:,:] = -gaussian2D(x,y,Lx//2,Ly//2)*corr + gaussian2D(x,y,Lx//6,Ly//6)*corr

########################################################
#COMPUTATION

if integrator == 'Euler':
    print('Integration : Euler')
    if rotating_frame == False:
        f = 0
        print('Coriolis parameter : ',f)
        u1,v1,eta1,u2,v2,eta2=get_SW_2layer_euler_R(t,x,y,Ny,u1,v1,eta1,u2,v2,eta2,detadt,dx,dy,dt,g,g2,H1,H2,f,obstacle)
    elif rotating_frame == True:
        print('Coriolis parameter : ',f)
        u1,v1,eta1,u2,v2,eta2=get_SW_2layer_euler_R(t,x,y,Ny,u1,v1,eta1,u2,v2,eta2,detadt,dx,dy,dt,g,g2,H1,H2,f,obstacle)

elif integrator == 'LP':
    print('Integration : Leap-Frog')
    if rotating_frame == False:
        f=0
        print('Coriolis parameter : ',f)
        u1,v1,eta1,u2,v2,eta2 = get_SW_2layer_LP_R(t,x,y,Ny,u1,v1,eta1,u2,v2,eta2,detadt,dx,dy,dt,g,g2,H1,H2,f,obstacle)
    elif rotating_frame == True:
        print('Coriolis parameter : ',f)
        u1,v1,eta1,u2,v2,eta2 = get_SW_2layer_LP_R(t,x,y,Ny,u1,v1,eta1,u2,v2,eta2,detadt,dx,dy,dt,g,g2,H1,H2,f,obstacle)

elif integrator == 'LP_m':
    print('Integration : Leap-Frog/Euler')
    if rotating_frame == False:
        f = 0
        u1,v1,eta1,u2,v2,eta2 = get_SW_2layers_LP_mixte_R(t,x,y,Ny,Nt,u1,v1,eta1,u2,v2,eta2,detadt,dx,dy,dt,g,g2,H1,H2,f,obstacle)
    elif rotating_frame == True:
        u1,v1,eta1,u2,v2,eta2 = get_SW_2layers_LP_mixte_R(t,x,y,Ny,Nt,u1,v1,eta1,u2,v2,eta2,detadt,dx,dy,dt,g,g2,H1,H2,f,obstacle)

###############################################################
#PLOT

if view_1D==True:
    plt.figure()
    n=25
    eps=1 
    for i in range(len(t)):
        plt.plot(x,eta1[i,:,n]+H1[:,n]+H2[:,n],label=r'$\eta_1$')
        plt.plot(x,eta2[i,:,n]+H2[:,n],label=r'$\eta_2$')
        plt.xlabel(r'$x$')
        plt.ylabel(r'$H$')
        plt.ylim(np.min(H2[:,n])-eps,np.max(H2[:,n]+H1[:,n])+eps)
        plt.title(str(i)+'/'+str(len(t)-1)+r' at $y=$'+str(n))
        plt.legend()
        plt.pause(0.01)
        plt.clf()
    plt.plot(x,eta1[i,:,n]+H1[:,n]+H2[:,n],label=r'$\eta_1$')
    plt.plot(x,eta2[i,:,n]+H2[:,n],label=r'$\eta_2$')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$H$')
    plt.title(str(i)+'/'+str(len(t)-1)+r' at $x=$'+str(n))
        
        
        
if view_2D==True:
    fig,(ax)=plt.subplots(1,2,figsize=(15,7))
    for i in range(len(t)):
        fig.suptitle(str(i)+'/'+str(len(t)-1))

        fig1=ax[0].pcolormesh(x,y,eta1[i,:,:],cmap=cmap,vmin=np.min(eta1),vmax=np.max(eta1))
        fig2=ax[1].pcolormesh(x,y,eta2[i,:,:],cmap=cmap,vmin=np.min(eta2),vmax=np.max(eta2))

        ax[0].set_xlabel(r'$x$')
        ax[0].set_ylabel(r'$y$')
        ax[0].set_title(r'$\eta_1$')
        ax[1].set_xlabel(r'$x$')
        ax[1].set_ylabel(r'$y$')
        ax[1].set_title(r'$\eta_2$')

        plt.pause(0.01)
        ax[0].clear()
        ax[1].clear()

if view_3D == True:
    fig, ax = plt.subplots(1,1,figsize=(10,7),subplot_kw={"projection": "3d"})
    xx,yy=np.meshgrid(x,y)
    for i in range(len(t)):

        fig.suptitle(str(i)+'/'+str(len(t)-1))
        fig3=ax.plot_surface(xx,yy,eta1[i,:,:]+H1+H2,cmap=cmap)
        ax.plot_surface(xx,yy,eta2[i,:,:]+H2,cmap=cmap)

        ax.set_zlim(np.min(eta1),np.max(eta2)+np.max(H2+H1))

        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        ax.set_zlabel('$H$')

        plt.pause(0.01)
        ax.clear()
    ax.plot_surface(xx,yy,eta1[-1,:,:]+H1+H2,cmap=cmap)
    ax.plot_surface(xx,yy,eta2[i,:,:]+H2,cmap=cmap)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_zlabel('$H$')
    ax.set_zlim(np.min(eta1),np.max(eta2)+np.max(H2+H1))

if view_cour == True:
    fig,(ax)=plt.subplots(1,2,figsize=(15,7))
    for i in range(len(t)):
        fig.suptitle(str(i)+'/'+str(len(t)-1))

        fig1=ax[0].pcolormesh(x,y,eta1[i,:,:],cmap=cmap,vmin=np.min(eta1),vmax=np.max(eta1))
        ax[0].quiver(x,y,u1[i,:,:],v1[i,:,:])
        fig2=ax[1].pcolormesh(x,y,eta2[i,:,:],cmap=cmap,vmin=np.min(eta2),vmax=np.max(eta2))
        ax[1].quiver(x,y,u2[i,:,:],v2[i,:,:])

        ax[0].set_xlabel(r'$x$')
        ax[0].set_ylabel(r'$y$')
        ax[0].set_title(r'$\eta_1$')
        ax[1].set_xlabel(r'$x$')
        ax[1].set_ylabel(r'$y$')
        ax[1].set_title(r'$\eta_2$')

        #fig.colorbar(fig1,ax=ax[0],label=r'$\eta_1$')

        plt.pause(0.01)
        ax[0].clear()
        ax[1].clear()
        
plt.show()


def save_nc(u1,v1,u2,v2,eta1,eta2):
    from netCDF4 import Dataset
    # Specify the netCDF file name
    netcdf_filename = 'output.nc'
    # Create a netCDF file in write mode
    with Dataset(netcdf_filename, 'w', format='NETCDF4') as ncfile:
        # Create dimensions matching the shape of your arrays ## common for all variables in t,x,y
        dim11 = ncfile.createDimension('time', eta1.shape[0])
        dim12 = ncfile.createDimension('x', eta1.shape[1])
        dim13 = ncfile.createDimension('y', eta1.shape[2])
        # Create variables and store the arrays
        var1 = ncfile.createVariable('eta1', 'f8', ('time','x', 'y'))
        var1[:] = eta1
        var2 = ncfile.createVariable('u1', 'f8', ('time','x', 'y'))
        var2[:] = u1  
        var3 = ncfile.createVariable('v1', 'f8', ('time','x', 'y'))
        var3[:] = v1    

        var11 = ncfile.createVariable('eta2', 'f8', ('time','x', 'y'))
        var11[:] = eta2
        var22 = ncfile.createVariable('u2', 'f8', ('time','x', 'y'))
        var22[:] = u2 
        var33 = ncfile.createVariable('v2', 'f8', ('time','x', 'y'))
        var33[:] = v2       

if save_netcf == True:
    save_nc(u1,v1,u2,v2,eta1,eta2)
