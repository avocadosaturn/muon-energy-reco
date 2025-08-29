import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D

ranges=[] #cm

f = h5py.File('/sdf/home/s/seohyeon/summer25/data/muon/mubar_1000_0-1GeV_gampixsim_5mmpp.h5')
pixels = f['pixel_hits']
coarse = f['coarse_hits']
meta = f['meta']



ev_mask_pixels = pixels['event id'] == 2
ev_pixels = pixels[ev_mask_pixels]

coords = []
for pixel in ev_pixels:
    hit_coord = np.array([pixel['pixel x'], pixel['pixel y'], pixel['hit z']])
    coords.append(hit_coord)


min_hits = 80
if len(coords) < min_hits:
    print(len(coords))
    print('not enough hits in this event')
    sys.exit()
else:
    pca = PCA(n_components = 3)
    pca_on_hits = pca.fit(coords)
    principal=pca.components_[0] #this is a unit vector


#getting projection of each coordinate onto the principal axis
parallel_proj={} #key is magnitude of projection, value is coordinate associated with it
antiparallel_proj = {}
for coord in coords:
    proj = np.dot(principal, coord) * principal
    mag_proj = np.linalg.norm(proj)

    cos = np.dot(principal, proj) / mag_proj #cos theta between the projection and principal axis. Will be 1 if parallel, -1 if antiparallel
    if cos > 0:
        parallel_proj[mag_proj] = coord
    elif cos < 0:
        antiparallel_proj[mag_proj] = coord

if len(parallel_proj)==0:
    track_length=max(antiparallel_proj)
    start=antiparallel_proj[min(antiparallel_proj)]
    stop=antiparallel_proj[max(antiparallel_proj)]
elif len(antiparallel_proj)==0:
    track_length=max(parallel_proj)
    start=parallel_proj[min(parallel_proj)]
    stop=parallel_proj[max(parallel_proj)]
else:
    track_length = max(parallel_proj)+max(antiparallel_proj)
    start=antiparallel_proj[max(antiparallel_proj)]
    stop=parallel_proj[max(parallel_proj)]


coords = np.vstack(coords)
coords = coords[~np.all(coords == start, axis=1)]
coords = coords[~np.all(coords == stop, axis=1)]

length_vector = track_length * principal
print(length_vector)



from matplotlib.widgets import Slider

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
origin=[0, 0, 0]

# plot axes
ax.quiver(*origin, 500, 0, 0, arrow_length_ratio = 0.01, color='black')
ax.text(500, 0, 0, 'x')
ax.quiver(*origin, 0, 500, 0, arrow_length_ratio=0.01, color='black')
ax.text(0, 500, 0, 'y')
ax.quiver(*origin, 0, 0, 500, arrow_length_ratio = 0.01, color='black')
ax.text(0, 0, 500, 'z')

ax.set_axis_off()

ax.set_xlim([-700, 700])
ax.set_ylim([-700, 700])
ax.set_zlim([-50, 1000])

        
# Plot points
ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], color='b', marker='o')
ax.scatter(*start, color='r', marker='o')
ax.scatter(*stop, color='r', marker='o')
# ax.quiver(0, 0, 0, *principal, length=800.0, arrow_length_ratio=0.01) #principal
ax.quiver(*start, *length_vector, arrow_length_ratio=0.01, color='r', linewidths=2)



plt.subplots_adjust(left=0.25, bottom=0.25)

#Slider axes: [left, bottom, width, height]
ax_elev = plt.axes([0.25, 0.1, 0.65, 0.03])
ax_azim = plt.axes([0.25, 0.15, 0.65, 0.03])

#Create sliders
slider_elev = Slider(ax_elev, 'Elevation', -90, 90, valinit=60)
slider_azim = Slider(ax_azim, 'Azimuth', -180, 180, valinit=-60)

#Update function to change view
def update(val):
   ax.view_init(elev=slider_elev.val, azim=slider_azim.val)
   fig.canvas.draw_idle()

slider_elev.on_changed(update)
slider_azim.on_changed(update)

plt.show()
plt.savefig(f'Images/muon_event{ev_n}_pixels.png')


