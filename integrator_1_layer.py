#######################################################################
# EULER                                                               #
#######################################################################
#Pour le Euler il faut garder un nombre Nt = 250 !
#dt = 6 # dtmax euler
#time = 1500 #max euler

def get_SW_euler_R(t,x,y,Ny,u,v,eta,dx,dy,dt,g,H,f,obstacle):
    for k in range(len(t)-1):
        for l in range(len(x)-1):
            for j in range(len(y)-1):
                u[k,0,:] = 0
                u[k,-1,:] = 0
                v[k,:,0] = 0
                v[k,:,-1] = 0

                if obstacle == 'detroit':
                    ####### COTE GAUCHE
                    u[k,Ny//3+5:Ny//3+15,0:24] = 0
                    v[k,Ny//3+5:Ny//3+15,0:24] = 0

                    ####### COTE DROITE
                    u[k,Ny//3+5:Ny//3+15,26:50] = 0
                    v[k,Ny//3+5:Ny//3+15,26:50] = 0

                elif obstacle == 'ile':
                    u[k,Ny//3+5:Ny//3+15,20:30] = 0
                    v[k,Ny//3+5:Ny//3+15,20:30] = 0
    
                u[k+1,l,j] = u[k,l,j]-dt*(g* (eta[k,l+1,j]-eta[k,l-1,j])/(2*dx) + f*v[k,l,j])
                v[k+1,l,j] = v[k,l,j]-dt*(g* (eta[k,l,j+1]-eta[k,l,j-1])/(2*dy) - f*u[k,l,j])
                eta[k+1,l,j] = eta[k,l,j]-dt*(H[l,j]*(u[k,l+1,j]-u[k,l-1,j])/(2*dx) + (v[k,l,j+1]-v[k,l,j-1])/(2*dy))        
    
    return u,v,eta

#######################################################################
# LEAP-FROG                                                           #
#######################################################################
#Pour le LP il faut garder un nombre Nt = 250 !
#dt = 20 #dt max LP
#time = 5000 #max LP

def get_SW_LP_R(t,x,y,Ny,u,v,eta,dx,dy,dt,g,H,f,obstacle):
    tspawn = 2
    for k in range(len(t)-1):
        for l in range(len(x)-1):
            for j in range(len(y)-1):
                u[k,0,:] = 0
                u[k,-1,:] = 0
                v[k,:,0] = 0
                v[k,:,-1] = 0

                if obstacle == 'detroit':
                    ####### COTE GAUCHE
                    u[k,Ny//3+5:Ny//3+15,0:24] = 0
                    v[k,Ny//3+5:Ny//3+15,0:24] = 0

                    ####### COTE DROITE
                    u[k,Ny//3+5:Ny//3+15,26:50] = 0
                    v[k,Ny//3+5:Ny//3+15,26:50] = 0

                elif obstacle == 'ile':
                    u[k,Ny//3+5:Ny//3+15,20:30] = 0
                    v[k,Ny//3+5:Ny//3+15,20:30] = 0

                #####################################################
                #SCHEMA MIXTE : 1) Euler ; 2) LP ; 3) Euler ; 4) LP.
                if t[k] <= tspawn:
                    #print('Euler init')
                    u[k+1,l,j] = u[k,l,j]-dt*(g* (eta[k,l+1,j]-eta[k,l-1,j])/(2*dx) + f*v[k,l,j])
                    v[k+1,l,j] = v[k,l,j]-dt*(g* (eta[k,l,j+1]-eta[k,l,j-1])/(2*dy) - f*u[k,l,j])
                    eta[k+1,l,j] = eta[k,l,j]-dt*(H[l,j]*((u[k,l+1,j]-u[k,l-1,j])/(2*dx) + (v[k,l,j+1]-v[k,l,j-1])/(2*dy)))

                elif (t[k] > tspawn):
                    #print('Leap-Frog')
                    u[k+1,l,j] = u[k-1,l,j]-2*dt*(g* (eta[k,l+1,j]-eta[k,l-1,j])/(2*dx) + f*v[k,l,j])
                    v[k+1,l,j] = v[k-1,l,j]-2*dt*(g* (eta[k,l,j+1]-eta[k,l,j-1])/(2*dy) - f*u[k,l,j])
                    eta[k+1,l,j] = eta[k-1,l,j]-2*dt*(H[l,j]*((u[k,l+1,j]-u[k,l-1,j])/(2*dx) + (v[k,l,j+1]-v[k,l,j-1])/(2*dy)))
    
    return u,v,eta


#######################################################################
# LEAP-FROG MIXTE EULER                                               #
#######################################################################
#Pour le LP mixte il faut garder un nombre Nt = 250 !
#dt = 20 #dt max LP_m
#time = 5000 #max LP_m

def get_SW_LP_mixte_R(t,x,y,Ny,Nt,u,v,eta,dx,dy,dt,g,H,f,obstacle):
    tspawn = 2
    tmixt = Nt//2
    dmixt = 50 #20
    for k in range(len(t)-1):
        for l in range(len(x)-1):
            for j in range(len(y)-1):
                u[k,0,:] = 0
                u[k,-1,:] = 0
                v[k,:,0] = 0
                v[k,:,-1] = 0

                if obstacle == 'detroit':
                    ####### COTE GAUCHE
                    u[k,Ny//3+5:Ny//3+15,0:24] = 0
                    v[k,Ny//3+5:Ny//3+15,0:24] = 0

                    ####### COTE DROITE
                    u[k,Ny//3+5:Ny//3+15,26:50] = 0
                    v[k,Ny//3+5:Ny//3+15,26:50] = 0

                elif obstacle == 'ile':
                    u[k,Ny//3+5:Ny//3+15,20:30] = 0
                    v[k,Ny//3+5:Ny//3+15,20:30] = 0

                #####################################################
                #SCHEMA MIXTE : 1) Euler ; 2) LP ; 3) Euler ; 4) LP.
                if t[k] <= tspawn:
                    #print('Euler init')
                    u[k+1,l,j] = u[k,l,j]-dt*(g* (eta[k,l+1,j]-eta[k,l-1,j])/(2*dx) + f*v[k,l,j])
                    v[k+1,l,j] = v[k,l,j]-dt*(g* (eta[k,l,j+1]-eta[k,l,j-1])/(2*dy) - f*u[k,l,j])
                    eta[k+1,l,j] = eta[k,l,j]-dt*(H[l,j]*((u[k,l+1,j]-u[k,l-1,j])/(2*dx) + (v[k,l,j+1]-v[k,l,j-1])/(2*dy)))

                elif (t[k] > tspawn) & (t[k]<=tmixt):
                    #print('Leap-Frog')
                    u[k+1,l,j] = u[k-1,l,j]-2*dt*(g* (eta[k,l+1,j]-eta[k,l-1,j])/(2*dx) + f*v[k,l,j])
                    v[k+1,l,j] = v[k-1,l,j]-2*dt*(g* (eta[k,l,j+1]-eta[k,l,j-1])/(2*dy) - f*u[k,l,j])
                    eta[k+1,l,j] = eta[k-1,l,j]-2*dt*(H[l,j]*((u[k,l+1,j]-u[k,l-1,j])/(2*dx) + (v[k,l,j+1]-v[k,l,j-1])/(2*dy)))
                
                elif (t[k] > tmixt) & (t[k]<=tmixt+dmixt):
                    #print('Euler')
                    u[k+1,l,j] = u[k,l,j]-dt*(g* (eta[k,l+1,j]-eta[k,l-1,j])/(2*dx) + f*v[k,l,j])
                    v[k+1,l,j] = v[k,l,j]-dt*(g* (eta[k,l,j+1]-eta[k,l,j-1])/(2*dy) - f*u[k,l,j])
                    eta[k+1,l,j] = eta[k,l,j]-dt*(H[l,j]*((u[k,l+1,j]-u[k,l-1,j])/(2*dx) + (v[k,l,j+1]-v[k,l,j-1])/(2*dy)))
                
                elif t[k] > tmixt+dmixt:
                    #print('Leap-Frog')
                    u[k+1,l,j] = u[k-1,l,j]-2*dt*(g* (eta[k,l+1,j]-eta[k,l-1,j])/(2*dx) + f*v[k,l,j])
                    v[k+1,l,j] = v[k-1,l,j]-2*dt*(g* (eta[k,l,j+1]-eta[k,l,j-1])/(2*dy) - f*u[k,l,j])
                    eta[k+1,l,j] = eta[k-1,l,j]-2*dt*(H[l,j]*((u[k,l+1,j]-u[k,l-1,j])/(2*dx) + (v[k,l,j+1]-v[k,l,j-1])/(2*dy)))
                
    return u,v,eta
