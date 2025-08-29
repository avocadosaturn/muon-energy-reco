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
    errors2=[]
    bin_indices = np.digitize(xpts, bins)  #finds corresponding bin for each xpt
    
    for i in range(1, len(bins)):
        bin_data = ypts[bin_indices == i]
        if len(bin_data) == 0:
            means.append(0)
            errors.append(0)
        else:
            means.append(bin_data.mean())
            errors.append(bin_data.std())


    rerrors = np.array(errors) / np.array(means)

    return bin_midpts, means, errors, rerrors
