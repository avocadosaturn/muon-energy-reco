import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.interpolate import splprep, splev
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D


f = h5py.File('/sdf/home/s/seohyeon/summer25/test_data/muon/muon1000ev_0-1gev_gampix_5mmpp.h5')
pixels = f['pixel_hits']
coarse = f['coarse_hits']
meta = f['meta']

ev_n = 2

ev_mask_pixels = pixels['event id'] == ev_n
ev_pixels = pixels[ev_mask_pixels]

points = np.array([[hit['pixel x'], hit['pixel y'], hit['hit z']] for hit in ev_pixels])

# tck, u = splprep([points[:, 0], points[:, 1], points[:, 2]], s=25150000)  # s is smoothing factor; s=0 gives interpolation
spl, u = make_splprep([points[:, 0], points[:, 1], points[:, 2]], k=3)

# u_fine = np.linspace(0, 1, 5000)
u = np.linspace(bspline.t[bspline.k], bspline.t[-bspline.k-1], 5000)
curve = bspline(u)
# x_fine, y_fine, z_fine = splev(u_fine, tck)


fig = plt.figure(figsize=(12, 7))
ax = fig.add_subplot(111, projection='3d')


ax.scatter(points[:,0], points[:,1], points[:,2], 'o', s=1, alpha =0.05, color='lightgray')
ax.plot(curve[:, 0], curve[:, 1], curve[:, 2], y_fine, z_fine, 'b-', )

# sliders
#Slider axes: [left, bottom, width, height]
ax_elev = plt.axes([0.25, 0.1, 0.65, 0.03])
ax_azim = plt.axes([0.25, 0.15, 0.65, 0.03])

#Create sliders
slider_elev = Slider(ax_elev, 'Elevation', -90, 90, valinit=24.7)
slider_azim = Slider(ax_azim, 'Azimuth', -180, 180, valinit=155)

#Update function to change view
def update(val):
   ax.view_init(elev=slider_elev.val, azim=slider_azim.val)
   fig.canvas.draw_idle()

slider_elev.on_changed(update)
slider_azim.on_changed(update)

plt.show()