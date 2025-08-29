def bin_data(xpts, ypts, n_bins):
    import numpy as np
    import math
    import sys

    bins = np.linspace(min(xpts), max(xpts), n_bins + 1) #making bins
    
    # calculating midpoints of bins
    bin_midpts = []
    for i in range(len(bins) - 1):
        x = (bins[i] + bins[i + 1]) / 2
        bin_midpts.append(x)
    
    
    # binning data
    means=[]
    errors=[]
    bin_indices = np.digitize(xpts, bins)  #finds corresponding bin for each xpt
    
    for i in range(1, len(bins)):
        bin_data = ypts[bin_indices == i]
        if len(bin_data) == 0:
            means.append(0)
            errors.append(0)
        else:
            means.append(bin_data.mean())
            errors.append(bin_data.std())
    
    
    rerrors = np.array(errors) / (means)



    # bins = np.linspace(min(xpts), max(xpts), n_bins + 1) #making bins
    
    # # calculating midpoints of bins
    # bin_midpts = []
    # for i in range(len(bins) - 1):
    #     x = (bins[i] + bins[i + 1]) / 2
    #     bin_midpts.append(x)

    # # binning data
    # binned_data = {i: [] for i in range(len(bins) - 1)} #initialize binning library
    # bin_indices = np.digitize(xpts, bins) - 1 #finds corresponding bin for each xpt
    # for i, bin_index in enumerate(bin_indices):
    #     if bin_index < 0:
    #         continue
    #     else:
    #         binned_data[bin_index].append(ypts[i]) #bins ypt into correct bin

    # # calculating averages of each bin
    # avgs = []
    # for i in binned_data:
    #     y_list = binned_data[i]
    #     avg_num = 0
    #     n_pts = len(y_list)

    #     if n_pts == 0:
    #         avg=0
    #         print('no points in this bin')
    #         avgs.append(avg)
    #     else:
    #         for y in y_list:
    #             avg_num += y
    #         avg = avg_num / n_pts
    #         avgs.append(avg)


    # # calculate RMS
    # errors = []
    # for i in binned_data:
    #     y_list = binned_data[i]
    #     error_num = 0
    #     n_pts = len(y_list)
        
    #     for y in y_list:
    #         error_num += abs(y - avgs[i]) ** 2

    #     if n_pts == 0:
    #         error=0
    #         errors.append(error)
    #     else:
    #         error = math.sqrt(error_num / n_pts)
    #         errors.append(error)


    # # calculate relative RMS
    # relative_rms = []
    # for i in range(len(bins) - 1):
    #     if avgs[i] == 0:
    #         relative_rms.append(0)
    #     else:
    #         relative_rms.append(errors[i] / avgs[i])

    return bin_midpts, means, errors, rerrors
