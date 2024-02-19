##################################
# SHALLOW WATER MODEL : 2 LAYERS #
#         Dimitri MOREAU         #
##################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#OPTIONAL : PRESENTATION

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['font.size'] = 15

##################################
#GENERAL PARAMETERS

omega = 2*np.pi/(24*3600)
f = 2*omega*np.sin(45)
#f=1e-2 #pour tester
g = 9.81

rho2 = 1020
rho1 = 1010
rho0 = 1013
g2= g*(rho1 - rho2)/rho0 #reduced gravity
print('Reduced gravity : ',g2)

H1=1
H2=1

rotating_frame=True
tide=False

cmap='RdBu_r'

###################################
#PLOT SELECTION

view_1D = False
view_2D = False
view_3D = True

##################################
#CHOICE OF PERTURBATION

IC='eta_detroit' #'decent', 'bidecent', 'cent','quadri'
amp_pert=0.05 #amplitude of the perturbation
detroit = False

###############################
#SPATIAL AND TIME DISCRETISATION
#DECLARATION OF VARIABLES

Lx = 5*1e3
Ly = 5*1e3

Nx = 50
Ny = 50
dx = Lx/Nx
dy = Ly/Ny

dt = 6
time = 1500 #max

Nt = time/dt
print('Nombre de cellules temporelles : ',Nt)

x=np.arange(0,Lx,dx)
y=np.arange(0,Ly,dy)
t=np.arange(0,time,dt)

u1=np.zeros([len(t),len(x),len(y)])
v1=np.zeros([len(t),len(x),len(y)])
eta1=np.zeros([len(t),len(x),len(y)])

u2=np.zeros([len(t),len(x),len(y)])
v2=np.zeros([len(t),len(x),len(y)])
eta2=np.zeros([len(t),len(x),len(y)])

detadt=np.zeros([len(t),len(x),len(y)])

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

if tide==True:
    T_tide=10
    for i in range(len(t)-1):
        eta1[i,:,:]=1e-2*np.cos(t[i]*(2*np.pi)/T_tide)

corr = 1e6 #correction pour la gaussienne 

if IC == 'decent':
    eta1[0,:,:] = gaussian2D(x,y,Lx//6,Ly//6)*corr
elif IC == 'cent':
    eta1[0,:,:] = gaussian2D(x,y,Lx//2,Ly//2)*corr
elif IC == 'eta_detroit':
    eta1[0,:,:] = gaussian2D(x,y,Lx//2,Ly//6)*corr
elif IC == 'bidecent':
    eta1[0,:,:] = gaussian2D(x,y,Lx//2,Ly//2)*corr + gaussian2D(x,y,Lx//6,Ly//6)*corr

########################################################
#COMPUTATION

def get_SW_2layer_euler_NR(u1,v1,eta1,u2,v2,eta2,dx,dy,dt,g,g2,H1,H2,detroit):
    for k in range(len(t)-1):
        for l in range(len(x)-1):
            for j in range(len(y)-1):
                u1[k,0,:]=0
                u1[k,-1,:]=0
                v1[k,:,0]=0
                v1[k,:,-1]=0

                u2[k,0,:]=0
                u2[k,-1,:]=0
                v2[k,:,0]=0
                v2[k,:,-1]=0

                if detroit==True:
                    ####### COTE GAUCHE
                    u1[k,Ny//3+5:Ny//3+15,0:24] = 0
                    v1[k,Ny//3+5:Ny//3+15,0:24] = 0
                    u2[k,Ny//3+5:Ny//3+15,0:24] = 0
                    v2[k,Ny//3+5:Ny//3+15,0:24] = 0

                    ####### COTE DROITE
                    u1[k,Ny//3+5:Ny//3+15,26:50] = 0
                    v1[k,Ny//3+5:Ny//3+15,26:50] = 0
                    u2[k,Ny//3+5:Ny//3+15,26:50] = 0
                    v2[k,Ny//3+5:Ny//3+15,26:50] = 0
    
                    
                    #1st layer
                u1[k+1,l,j]=u1[k,l,j]-dt*(g* (eta1[k,l+1,j]-eta1[k,l-1,j])/(2*dx))
                v1[k+1,l,j]=v1[k,l,j]-dt*(g* (eta1[k,l,j+1]-eta1[k,l,j-1])/(2*dy))
                    
                detadt[k,l,j]=-H2*((u2[k,l+1,j] - u2[k,l-1,j])/(2*dx) + (v2[k,l,j+1] - v2[k,l,j-1])/(2*dy))
                    
                eta1[k+1,l,j]=eta1[k,l,j]-dt*(H1*(u1[k,l+1,j]-u1[k,l-1,j])/(2*dx) + (v1[k,l,j+1]-v1[k,l,j-1])/(2*dy)) - detadt[k,l,j]
                    
                    #2nd layer
                u2[k+1,l,j]=u2[k,l,j]-dt*(g* (eta1[k,l+1,j]-eta1[k,l-1,j])/(2*dx) - g2*(eta2[k,l+1,j]-eta2[k,l-1,j])/(2*dx))
                v2[k+1,l,j]=v2[k,l,j]-dt*(g* (eta1[k,l,j+1]-eta1[k,l,j+1])/(2*dy) - g2*(eta2[k,l,j+1]-eta2[k,l,j+1])/(2*dy))
                eta2[k+1,l,j]=eta2[k,l,j]-dt*H2*((u2[k,l+1,j]-u2[k,l-1,j])/(2*dx) + (v2[k,l,j+1]-v2[k,l,j-1])/(2*dy))
                        
    return u1,v1,eta1,u2,v2,eta2

def get_SW_2layer_euler_R(u1,v1,eta1,u2,v2,eta2,dx,dy,dt,g,g2,H1,H2,f,detroit):
    for k in range(len(t)-1):
        for l in range(len(x)-1):
            for j in range(len(y)-1):
                u1[k,0,:]=0
                u1[k,-1,:]=0
                v1[k,:,0]=0
                v1[k,:,-1]=0

                u2[k,0,:]=0
                u2[k,-1,:]=0
                v2[k,:,0]=0
                v2[k,:,-1]=0

                if detroit==True:
                    ####### COTE GAUCHE
                    u1[k,Ny//3+5:Ny//3+15,0:24] = 0
                    v1[k,Ny//3+5:Ny//3+15,0:24] = 0
                    u2[k,Ny//3+5:Ny//3+15,0:24] = 0
                    v2[k,Ny//3+5:Ny//3+15,0:24] = 0

                    ####### COTE DROITE
                    u1[k,Ny//3+5:Ny//3+15,26:50] = 0
                    v1[k,Ny//3+5:Ny//3+15,26:50] = 0
                    u2[k,Ny//3+5:Ny//3+15,26:50] = 0
                    v2[k,Ny//3+5:Ny//3+15,26:50] = 0
    
                
                #1st layer
                u1[k+1,l,j]=u1[k,l,j]-dt*(g* (eta1[k,l+1,j]-eta1[k,l-1,j])/(2*dx) + f*v1[k,l,j])
                v1[k+1,l,j]=v1[k,l,j]-dt*(g* (eta1[k,l,j+1]-eta1[k,l,j-1])/(2*dy) - f*u1[k,l,j])
                    
                detadt[k,l,j]=-H2*((u2[k,l+1,j] - u2[k,l-1,j])/(2*dx) + (v2[k,l,j+1] - v2[k,l,j-1])/(2*dy))
                    
                eta1[k+1,l,j]=eta1[k,l,j]-dt*(H1*(u1[k,l+1,j]-u1[k,l-1,j])/(2*dx) + (v1[k,l,j+1]-v1[k,l,j-1])/(2*dy)) - detadt[k,l,j]
                    
                #2nd layer
                u2[k+1,l,j]=u2[k,l,j]-dt*(g* (eta1[k,l+1,j]-eta1[k,l-1,j])/(2*dx) - g2*(eta2[k,l+1,j]-eta2[k,l-1,j])/(2*dx)+f*v2[k,l,j])
                v2[k+1,l,j]=v2[k,l,j]-dt*(g* (eta1[k,l,j+1]-eta1[k,l,j+1])/(2*dy) - g2*(eta2[k,l,j+1]-eta2[k,l,j+1])/(2*dy)-f*u2[k,l,j])
                eta2[k+1,l,j]=eta2[k,l,j]-dt*H2*((u2[k,l+1,j]-u2[k,l-1,j])/(2*dx) + (v2[k,l,j+1]-v2[k,l,j-1])/(2*dy))
                        
    return u1,v1,eta1,u2,v2,eta2

if rotating_frame==False:
    u1,v1,eta1,u2,v2,eta2=get_SW_2layer_euler_NR(u1,v1,eta1,u2,v2,eta2,dx,dy,dt,g,g2,H1,H2,detroit)
elif rotating_frame==True:
    u1,v1,eta1,u2,v2,eta2=get_SW_2layer_euler_R(u1,v1,eta1,u2,v2,eta2,dx,dy,dt,g,g2,H1,H2,f,detroit)

###############################################################
#PLOT

if view_1D==True:
    plt.figure()
    n=25
    eps=1
    
    for i in range(len(t)):
        plt.plot(x,eta1[i,:,n]+H1+H2,label=r'$\eta_1$')
        plt.plot(x,eta2[i,:,n]+H2,label=r'$\eta_2$')
        plt.xlabel(r'$x$')
        plt.ylabel(r'$H$')
        plt.ylim(H2-eps,H2+H1+eps)
        plt.title(str(i)+'/'+str(len(t))+r' at $x=$'+str(n))
        plt.legend()
        plt.pause(0.01)
        plt.clf()
        
        
        
if view_2D==True:
    fig,(ax)=plt.subplots(1,2,figsize=(15,7))
    for i in range(len(t)-1):
        fig.suptitle(str(i)+'/'+str(len(t)-1))

        fig1=ax[0].contourf(x,y,eta1[i,:,:],cmap=cmap,vmin=np.min(eta1),vmax=np.max(eta1))
        fig2=ax[1].contourf(x,y,eta2[i,:,:],cmap=cmap,vmin=np.min(eta2),vmax=np.max(eta2))

        plt.pause(0.01)
        ax[0].clear()
        ax[1].clear()


if view_3D == True:
    fig, ax = plt.subplots(1,1,figsize=(10,7),subplot_kw={"projection": "3d"})
    xx,yy=np.meshgrid(x,y)
    for i in range(len(t)-1):

        fig.suptitle(str(i)+'/'+str(len(t)-1))
        fig3=ax.plot_surface(xx,yy,eta1[i,:,:]+H1+H2,cmap=cmap)
        ax.plot_surface(xx,yy,eta2[i,:,:]+H2,cmap=cmap)

        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        ax.set_zlabel('$H$')

        plt.pause(0.01)
        ax.clear()
    ax.plot_surface(xx,yy,eta1[-1,:,:]+H1+H2,cmap=cmap)
    ax.plot_surface(xx,yy,eta2[i,:,:]+H2,cmap=cmap)
        
plt.show()






