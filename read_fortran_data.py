import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


Nt, Nx, Ny = 250, 50, 50

A = np.loadtxt('output.txt',skiprows=2)
eta = A.reshape([Nt, Nx, Ny], order='F')

zmin = float(np.min(eta))
zmax = float(np.max(eta))
zoff = zmin  # contour offset plane

x, y = np.linspace(0,1,Nx), np.linspace(0,1,Ny)

xx, yy = np.meshgrid(x,y)



fig, ax = plt.subplots(1,1,figsize=(10,7),subplot_kw={"projection": "3d"})
xx,yy=np.meshgrid(x,y)
for i in range(Nt):

	fig.suptitle(str(i)+'/'+str(Nt))

	fig3 = ax.plot_wireframe(xx,yy,eta[i,:,:]+0.,colors='k', rstride=2, cstride=2, linewidth=1)
	ax.contourf(xx, yy, eta[i,:,:],
           zdir='z',
           offset=eta.min()-0.3,
           cmap='RdBu_r')


	ax.set_zlim(-0.5,0.5)

	ax.set_xlabel('$x$')
	ax.set_ylabel('$y$')
	ax.set_zlabel('$H$')


	plt.pause(0.01)
	ax.clear()


fig3 = ax.plot_wireframe(xx,yy,eta[-1,:,:],colors='k', rstride=2, cstride=2, linewidth=1)
ax.contourf(xx, yy, eta[-1,:,:], zdir='z', offset=eta.min()-0.3, cmap="RdBu_r")

ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$H$')
ax.set_zlim(-0.5,0.5)



plt.show()
