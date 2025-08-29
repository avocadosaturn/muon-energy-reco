def bin20(xpts, ypts):
    import numpy as np
    import math
    import sys

    bins = np.linspace(0, math.ceil(max(xpts)), 21) #making bins
    
    # calculating midpoints of bins
    bin_midpts = []
    for i in range(len(bins) - 1):
        x = (bins[i] + bins[i + 1]) / 2
        bin_midpts.append(x)

    # binning data
    binned_data = {i: [] for i in range(len(bins) - 1)} #initialize binning library
    bin_indices = np.digitize(xpts, bins) - 1 #finds corresponding bin for each xpt
    for i, bin_index in enumerate(bin_indices):
        binned_data[bin_index].append(ypts[i]) #bins ypt into correct bin

    # calculating averages of each bin
    avgs = []
    for i in binned_data:
        y_list = binned_data[i]
        avg_num = 0
        n_pts = len(y_list)

        if n_pts == 0:
            avg = 0
            print('no points in this bin')
            avgs.append(avg)
        else:
            for y in y_list:
                avg_num += y
            avg = avg_num / n_pts
            avgs.append(avg)


    # calculate RMS
    errors = []
    for i in binned_data:
        y_list = binned_data[i]
        error_num = 0
        n_pts = len(y_list)
        
        for y in y_list:
            error_num += abs(y - avgs[i]) ** 2

        error = math.sqrt(error_num / n_pts)
        errors.append(error)


    # calculate relative RMS
    relative_rms = []
    for i in range(len(bins) - 1):
        relative_rms.append(errors[i] / avgs[i])


    return bin_midpts, avgs, errors, relative_rms
