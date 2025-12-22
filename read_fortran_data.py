import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import matplotlib.font_manager as fm
font_path = "/home/dmoreau/Bureau/courbd.ttf"
prop = fm.FontProperties(fname=font_path)



Nt, Nx, Ny = 250, 50, 50

A = np.loadtxt('output.txt',skiprows=2)
eta = A.reshape([Nt, Nx, Ny], order='F')

zmin = float(np.min(eta))
zmax = float(np.max(eta))
zoff = zmin  # contour offset plane

x, y = np.linspace(0,5000,Nx), np.linspace(0,5000,Ny)

xx, yy = np.meshgrid(x,y)



fig, ax = plt.subplots(1,1,figsize=(10,7),subplot_kw={"projection": "3d"})

norm = mpl.colors.Normalize(vmin=-0.15, vmax=0.15)
cmap = plt.get_cmap("RdBu_r")

mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
mappable.set_array([])
cbar = fig.colorbar(mappable, ax=ax, shrink=0.8, pad=0.1, extend='both')
cbar.set_label('η', fontproperties=prop, fontsize=13)


# Colorbar ticks font
for tick in cbar.ax.get_yticklabels():
	tick.set_fontproperties(prop)
	tick.set_fontsize(11)
	tick.set_fontweight('bold')

#cbar.ax.tick_params(direction='in',length=4,width=1.)
#cbar.ax.yaxis.set_ticks_position('both')
#cbar.ax.yaxis.set_tick_params(labelleft=True, labelright=True)





for i in range(Nt):

	fig.suptitle('t = '+str(i)+'/'+str(Nt), fontproperties=prop, fontsize=15)

	fig3 = ax.plot_wireframe(xx,yy,eta[i,:,:]+0.,colors='k', rstride=2, cstride=2, linewidth=1)
	ax.contourf(xx, yy, eta[i,:,:],
           zdir='z',
           offset=eta.min()-0.3,
           cmap=cmap, norm=norm)


	ax.set_zlim(-0.5,0.5)

	ax.set_xlabel('x', fontproperties=prop, fontsize=11)
	ax.set_ylabel('y', fontproperties=prop, fontsize=11)
	ax.set_zlabel('η', fontproperties=prop, fontsize=11)


	for spine in ax.spines.values():
		spine.set_linewidth(2)
	for tick in ax.get_xticklabels():
		tick.set_fontweight('bold')
	for tick in ax.get_yticklabels():
		tick.set_fontweight('bold')
	for tick in ax.get_zticklabels():
		tick.set_fontweight('bold')

	for label in ax.get_xticklabels() + ax.get_yticklabels() + ax.get_zticklabels():
		label.set_fontproperties(prop)
		label.set_fontsize(11)




	plt.pause(0.01)
	ax.clear()


fig3 = ax.plot_wireframe(xx,yy,eta[-1,:,:],colors='k', rstride=2, cstride=2, linewidth=1)
ax.contourf(xx, yy, eta[-1,:,:], zdir='z', offset=eta.min()-0.3, cmap=cmap,norm=norm)

ax.set_xlabel('x', fontproperties=prop, fontsize=13)
ax.set_ylabel('y', fontproperties=prop, fontsize=13)
ax.set_zlabel('η', fontproperties=prop, fontsize=13)
ax.set_zlim(-0.5,0.5)



plt.show()
