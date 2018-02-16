#!/usr/bin/env python

import h5py
import numpy as np
import matplotlib
matplotlib.rcParams["backend"] = "TkAgg"
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx

n_freq_points = 250
filename = 'ring_resonator-freq-monitor.h5'
f = h5py.File(filename, 'r')

# Get the data
field = np.array(f[list(f.keys())[0]])

# Get monitor metadata
max_freq = f.attrs.get('maxfreq')
min_freq = f.attrs.get('minfreq')
dt = f.attrs.get('dt')
xmin = f.attrs.get('xmin')
xmax = f.attrs.get('xmax')
ymin = f.attrs.get('ymin')
ymax = f.attrs.get('ymax')
f.close()


# For weakly resonant structures sometimes you sometimes need to
# apodize the field to avoid having the source completely dominate
# the output.
# Note: Apodization completely messes up any energy normalization.
#
# Ring resonator examples needs no apodization
t0=0.   
tau=0.01 # small tau means a fast ramp up
t = np.arange(0, field.shape[-1])*dt
apodize = 1 / (1 + np.exp(-(t-t0)/tau))
field = field * apodize

# You may need to pad with zeros or cut the input fields to FFT 
# in order to get the desired amount of frequency points on the output. 
n_sample_points = 2 * int(n_freq_points/(max_freq-min_freq) * max_freq)

fft_field = abs(np.fft.fft(field, n=n_sample_points, axis=-1).real)
fft_field = fft_field[:, :, :fft_field.shape[-1]//2]
freqs = np.fft.fftfreq(n_sample_points, d=dt)[:(n_sample_points/2)]

# We are only interested on the frequencies inside the source bandwidth
mask = (freqs >= min_freq) * (freqs <= max_freq)
freqs = freqs[mask]
fft_field = fft_field[:, :, mask]


# Plot Results
central_freq = (max_freq+min_freq)*0.5

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, bottom=0.25)
field_plot = plt.imshow(fft_field[:,::-1, find_nearest(freqs, central_freq)].T,
	                    extent=[xmin, xmax, ymin, ymax])
plt.colorbar()

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.1, 0.1, 0.65, 0.03], facecolor=axcolor)
sfreq = Slider(axfreq, 'Freq', np.min(freqs), np.max(freqs), valinit=central_freq)

def update(val):
    freq = sfreq.val
    field_plot.set_data(fft_field[:,::-1, find_nearest(freqs, freq)].T)
    field_plot.autoscale()
    fig.canvas.draw_idle()
sfreq.on_changed(update)

plt.show()
