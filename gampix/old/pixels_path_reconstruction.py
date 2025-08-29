import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D


f = h5py.File('/sdf/home/s/seohyeon/summer25/data/muon/mubar_1000_0-1GeV_gampixsim_5mmpp.h5')
pixels = f['pixel_hits']
meta = f['meta']
ev_n = 856

ev_mask_pixels = pixels['event id'] == ev_n
ev_pixels = pixels[ev_mask_pixels]

#getting the coordinates of each hit
coords = np.array([[hit['pixel x'], hit['pixel y'], hit['hit z']] for hit in ev_pixels])
print(f'number of pts: {len(coords)}')

# pca analysis
min_hits = 40
if len(coords) < min_hits:
    print(f'number of hits: {len(coords)}')
    print('not enough hits in this event')
    sys.exit()
else:
    pca = PCA(n_components = 3)
    pca_on_hits = pca.fit(coords)
    principal=pca.components_[0] #this is a unit vector


#getting projection of each coordinate onto the principal axis
projections ={}
for coord in coords:
    mag_proj = np.dot(principal, coord) #signed projection scalar
    projections[mag_proj] = coord


# calculate range
track_length = max(projections) - min(projections)
start = projections[min(projections)]
stop = projections[max(projections)]
print(f'track length: {track_length}')
length_vector = (stop - start) * track_length / np.linalg.norm(stop-start)


#removing start and stop pts
coords = np.vstack(coords)
coords = coords[~np.all(coords == start, axis=1)]
coords = coords[~np.all(coords == stop, axis=1)]


#----------------------------------------------------------
from matplotlib.widgets import Slider

fig = plt.figure(figsize=(12, 7))
ax = fig.add_subplot(111, projection='3d')
ax.set_title(f'event {ev_n}')

origin=[0, 0, 0]
max_z = max(coords[:, 2])

## !!!!!! x here is y on plot, y->z, z->x !!!!!!!!!
# plot limits
ax.view_init(elev=24.7, azim=155)
ax.set_xlim([-300, 300])
ax.set_zlim([-300, 300])
ax.set_ylim([-(max_z+50), max_z+50])

# axes
ax.set_axis_off()
ax.quiver(*origin, 300, 0, 0, arrow_length_ratio = 0.01, color='silver', alpha=0.5)
ax.text(300, 0, 0, 'y')
ax.quiver(*origin, 0, max_z + 50, 0, arrow_length_ratio=0.01, color='silver', alpha=0.5)
ax.text(0, max_z+50, 0, 'z')
ax.quiver(*origin, 0, 0, 300, arrow_length_ratio = 0.01, color='silver', alpha=0.5)
ax.text(0, 0, 300, 'x')

# anode
z = np.linspace(-300, 300, 10)
x = np.linspace(-300, 300, 10)
z, x = np.meshgrid(z, x)
y = np.zeros_like(z)
ax.plot_surface(x, y, z, color='lavender', alpha=0.1)

        
# plotting
ax.scatter(coords[:, 1], coords[:, 2], coords[:, 0], color='lightsteelblue', marker='o', alpha=0.1, s=10)
ax.scatter(start[1], start[2], start[0], color='r', marker='o', s=50, alpha=0.5)
ax.scatter(stop[1], stop[2], stop[0], color='r', marker='o', s=50, alpha=0.5)
# ax.quiver(0, 0, 0, *principal, length=800.0, arrow_length_ratio=0.01) #principal
ax.quiver(start[1], start[2], start[0], length_vector[1], length_vector[2], length_vector[0], arrow_length_ratio=0.01, color='blue', linewidths=2)


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


plt.savefig(f'Images/muon_event{ev_n}_pixels.png')
plt.show()



