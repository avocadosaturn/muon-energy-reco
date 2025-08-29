import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import binned_statistic_2d
from scipy.interpolate import splprep, splev


f = h5py.File('/sdf/home/s/seohyeon/summer25/test_data/muon/muon1000ev_0-1gev_gampix_5mmpp.h5')
pixels = f['pixel_hits']
coarse = f['coarse_hits']
meta = f['meta']

ev_n = 9

ev_mask_pixels = pixels['event id'] == ev_n
ev_pixels = pixels[ev_mask_pixels]


points = np.array([[hit['pixel x'], hit['pixel y'], hit['hit z']] for hit in ev_pixels])


bins1d = 5
statistic, x_edges, y_edges, binnumber = binned_statistic_2d(
    points[:, 0], points[:, 1], points[:, 2], statistic='mean', bins=(bins1d, bins1d)
)

KE = meta[ev_n]['primary energy']
print(f'number of points: {len(points)}')
print(f'energy: {KE}')


# pca ------------
min_hits = 80
if len(points) < min_hits:
    print(len(points))
    print('not enough hits in this event')
    sys.exit()
else:
    pca = PCA(n_components = 3)
    pca_on_hits = pca.fit(points)
    principal=pca.components_[0] #this is a unit vector


projections ={}
for pt in points:
    mag_proj = np.dot(principal, pt) #signed projection scalar
    projections[mag_proj] = pt


start = projections[min(projections)]
stop = projections[max(projections)]
start = np.array(start)
stop = np.array(stop)

naive_length = max(projections) - min(projections)
naive_vector = (stop - start) * naive_length / np.linalg.norm(stop - start)

# binning and plotting --------------
# settings
fig = plt.figure(figsize=(12, 7))
ax = fig.add_subplot(111, projection='3d')
ax.set_title(f'event {ev_n}')

ax.view_init(elev=24.7, azim=155)

# max_z = max(points[:, 2])
# ax.set_xlim([-300, 300])
# ax.set_ylim([-300, 300])
# ax.set_zlim([-(max_z+5), max_z+50])


#binning
spl_pts = []
for i, row in enumerate(statistic):
    for j, val in enumerate(row):
        if np.isnan(val)==True:
            continue
        else:
            
            i_prime = i + 1
            j_prime = j + 1
            bins1d_prime = bins1d + 2
            bin_n = (i_prime * bins1d_prime) + j_prime

            indices = np.where(binnumber == bin_n)

            x_avg = 0
            y_avg = 0
            n = 0
            for index in indices[0]:
                x_avg += points[index][0]
                y_avg += points[index][1]
                n = n+1
            x_avg = x_avg / n
            y_avg = y_avg / n

            spl_pt = [x_avg, y_avg, val]
            spl_pts.append(spl_pt)
            ax.scatter(x_avg, y_avg, val, color='black', s=6)

# sort binned points for spline
sort_spl_pts = {}
for pt in spl_pts:
    dist = np.linalg.norm(pt - start)
    sort_spl_pts[dist] = pt
spl_pts = []
for key in sorted(sort_spl_pts.keys()):
    spl_pts.append(sort_spl_pts[key])
spl_pts.insert(0, start)
spl_pts.append(stop)
spl_pts = np.array(spl_pts)

# spline
tck, u = splprep([spl_pts[:, 0], spl_pts[:, 1], spl_pts[:, 2]], s=50)
u_fine = np.linspace(0, 1, 1000)
x_fine, y_fine, z_fine = splev(u_fine, tck)


ax.scatter(points[:,0], points[:,1], points[:,2], 'o', s=20, color='lightpink', alpha = 0.01) # gampix points
plt.plot(x_fine, y_fine, z_fine, color='black') # spline curve
ax.scatter(*start, 'o', color='darkblue', s=10) 
ax.scatter(*stop, 'o', color='darkblue', s=10)
ax.quiver(*start, *naive_vector, arrow_length_ratio = 0.01, color='darkblue', linewidths=2)



# sliders
#slider axes: [left, bottom, width, height]
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











