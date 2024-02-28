#######################################################################
# EULER                                                               #
#######################################################################
#Pour le Euler il faut garder un nombre Nt = 250 !
#dt = 6 # dtmax euler
#time = 1500 #max euler

def get_SW_2layer_euler_R(t,x,y,Ny,u1,v1,eta1,u2,v2,eta2,detadt,dx,dy,dt,g,g2,H1,H2,f,obstacle):
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

                if obstacle == 'detroit':
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
                elif obstacle == 'ile':
                    u1[k,Ny//3+5:Ny//3+15,20:30] = 0
                    v1[k,Ny//3+5:Ny//3+15,20:30] = 0
                    u2[k,Ny//3+5:Ny//3+15,20:30] = 0
                    v2[k,Ny//3+5:Ny//3+15,20:30] = 0
    
                
                #1st layer
                u1[k+1,l,j]=u1[k,l,j]-dt*(g* (eta1[k,l+1,j]-eta1[k,l-1,j])/(2*dx) + f*v1[k,l,j])
                v1[k+1,l,j]=v1[k,l,j]-dt*(g* (eta1[k,l,j+1]-eta1[k,l,j-1])/(2*dy) - f*u1[k,l,j])
                detadt[k,l,j]=-H2[l,j]*((u2[k,l+1,j] - u2[k,l-1,j])/(2*dx) + (v2[k,l,j+1] - v2[k,l,j-1])/(2*dy))
                eta1[k+1,l,j]=eta1[k,l,j]-dt*(H1[l,j]*(u1[k,l+1,j]-u1[k,l-1,j])/(2*dx) + (v1[k,l,j+1]-v1[k,l,j-1])/(2*dy)) - detadt[k,l,j]
                    
                #2nd layer
                u2[k+1,l,j]=u2[k,l,j]-dt*(g* (eta1[k,l+1,j]-eta1[k,l-1,j])/(2*dx) - g2*(eta2[k,l+1,j]-eta2[k,l-1,j])/(2*dx)+f*v2[k,l,j])
                v2[k+1,l,j]=v2[k,l,j]-dt*(g* (eta1[k,l,j+1]-eta1[k,l,j+1])/(2*dy) - g2*(eta2[k,l,j+1]-eta2[k,l,j+1])/(2*dy)-f*u2[k,l,j])
                eta2[k+1,l,j]=eta2[k,l,j]-dt*H2[l,j]*((u2[k,l+1,j]-u2[k,l-1,j])/(2*dx) + (v2[k,l,j+1]-v2[k,l,j-1])/(2*dy))
                        
    return u1,v1,eta1,u2,v2,eta2

#######################################################################
# LEAP-FROG                                                           #
#######################################################################
#Pour le LP il faut garder un nombre Nt = 250 !
#dt = 20 #dt max LP
#time = 5000 #max LP

def get_SW_2layer_LP_R(t,x,y,Ny,u1,v1,eta1,u2,v2,eta2,detadt,dx,dy,dt,g,g2,H1,H2,f,obstacle):
    tspawn = 2
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

                if obstacle == 'detroit':
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
                elif obstacle == 'ile':
                    u1[k,Ny//3+5:Ny//3+15,20:30] = 0
                    v1[k,Ny//3+5:Ny//3+15,20:30] = 0
                    u2[k,Ny//3+5:Ny//3+15,20:30] = 0
                    v2[k,Ny//3+5:Ny//3+15,20:30] = 0
    
                if t[k] <= tspawn:
                    #1st layer
                    u1[k+1,l,j]=u1[k,l,j]-dt*(g* (eta1[k,l+1,j]-eta1[k,l-1,j])/(2*dx) + f*v1[k,l,j])
                    v1[k+1,l,j]=v1[k,l,j]-dt*(g* (eta1[k,l,j+1]-eta1[k,l,j-1])/(2*dy) - f*u1[k,l,j])
                    detadt[k,l,j]=-H2[l,j]*((u2[k,l+1,j] - u2[k,l-1,j])/(2*dx) + (v2[k,l,j+1] - v2[k,l,j-1])/(2*dy))
                    eta1[k+1,l,j]=eta1[k,l,j]-dt*(H1[l,j]*(u1[k,l+1,j]-u1[k,l-1,j])/(2*dx) + (v1[k,l,j+1]-v1[k,l,j-1])/(2*dy)) - detadt[k,l,j]
                        
                    #2nd layer
                    u2[k+1,l,j]=u2[k,l,j]-dt*(g* (eta1[k,l+1,j]-eta1[k,l-1,j])/(2*dx) - g2*(eta2[k,l+1,j]-eta2[k,l-1,j])/(2*dx)+f*v2[k,l,j])
                    v2[k+1,l,j]=v2[k,l,j]-dt*(g* (eta1[k,l,j+1]-eta1[k,l,j+1])/(2*dy) - g2*(eta2[k,l,j+1]-eta2[k,l,j+1])/(2*dy)-f*u2[k,l,j])
                    eta2[k+1,l,j]=eta2[k,l,j]-dt*H2[l,j]*((u2[k,l+1,j]-u2[k,l-1,j])/(2*dx) + (v2[k,l,j+1]-v2[k,l,j-1])/(2*dy))

                elif t[k] > tspawn:
                    #1st layer
                    u1[k+1,l,j]=u1[k-1,l,j]-2*dt*(g* (eta1[k,l+1,j]-eta1[k,l-1,j])/(2*dx) + f*v1[k,l,j])
                    v1[k+1,l,j]=v1[k-1,l,j]-2*dt*(g* (eta1[k,l,j+1]-eta1[k,l,j-1])/(2*dy) - f*u1[k,l,j])
                    detadt[k,l,j]=-H2[l,j]*((u2[k,l+1,j] - u2[k,l-1,j])/(2*dx) + (v2[k,l,j+1] - v2[k,l,j-1])/(2*dy))
                    eta1[k+1,l,j]=eta1[k-1,l,j]-2*dt*(H1[l,j]*(u1[k,l+1,j]-u1[k,l-1,j])/(2*dx) + (v1[k,l,j+1]-v1[k,l,j-1])/(2*dy)) - detadt[k,l,j]
                        
                    #2nd layer
                    u2[k+1,l,j]=u2[k-1,l,j]-2*dt*(g* (eta1[k,l+1,j]-eta1[k,l-1,j])/(2*dx) - g2*(eta2[k,l+1,j]-eta2[k,l-1,j])/(2*dx)+f*v2[k,l,j])
                    v2[k+1,l,j]=v2[k-1,l,j]-2*dt*(g* (eta1[k,l,j+1]-eta1[k,l,j+1])/(2*dy) - g2*(eta2[k,l,j+1]-eta2[k,l,j+1])/(2*dy)-f*u2[k,l,j])
                    eta2[k+1,l,j]=eta2[k-1,l,j]-2*dt*H2[l,j]*((u2[k,l+1,j]-u2[k,l-1,j])/(2*dx) + (v2[k,l,j+1]-v2[k,l,j-1])/(2*dy))


                            
    return u1,v1,eta1,u2,v2,eta2
