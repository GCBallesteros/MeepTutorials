#!/usr/bin/env python
import h5py
import numpy as np
import argparse

import matplotlib

matplotlib.rcParams["backend"] = "TkAgg"
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.colors import ListedColormap

# Parameters for epsilon
default_eps_parameters = {"interpolation": "spline36", "alpha": 1.0, "aspect": "equal"}


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx


def sigmoid_apodization(t, t0, tau):
    """Sigmoid apodization function.

    Used to attenuate the beginning of the signal.

    Arguments:
    t: Array with times.
    t0: Center point of the sigmoid.
    tau: Charectistic ramp up time. Bigger values produce
         smoother ramp ups."""
    return 1 / (1 + np.exp(-(t - t0) / tau))


def load_monitor(fname):
    f = h5py.File(fname, "r")

    field = np.array(f[list(f.keys())[0]])

    metadata = {
        "max_freq": f.attrs.get("maxfreq"),
        "min_freq": f.attrs.get("minfreq"),
        "dt": f.attrs.get("dt"),
        "xmin": f.attrs.get("xmin"),
        "xmax": f.attrs.get("xmax"),
        "ymin": f.attrs.get("ymin"),
        "ymax": f.attrs.get("ymax"),
        "res": f.attrs.get("res"),
    }

    f.close()

    return field, metadata


def load_struc(ename):
    # Load the epsilon
    f = h5py.File(ename, "r")
    eps = np.array(f[list(f.keys())[0]])
    f.close()
    eshape = np.array(eps.shape)

    # Shape the epsilon variable depending on dimensions
    if eps.ndim == 3:
        eps = eps[:, ::-1, int(eshape[2] / 2)].T
    else:
        eps = eps[:, ::-1].T

    return eps


def transform_field(field, metadata, npoints, apodization_func=None):
    max_freq = metadata["max_freq"]
    min_freq = metadata["min_freq"]
    dt = metadata["dt"]

    if apodization_func is not None:
        # For weakly resonant structures sometimes you sometimes need to
        # apodize the field to avoid having the source completely dominate
        # the output.
        # Note: Apodization  messes up energy normalization.
        t = np.arange(0, field.shape[-1]) * dt
        field = field * apodization_func(t)

    # Depending on how many frequencies you want to output you may need
    # to pad with zeros the input to the FFT.
    n_sample_points = 2 * int(npoints / (max_freq - min_freq) * max_freq)

    fft_field = abs(np.fft.fft(field, n=n_sample_points, axis=-1).real)
    fft_field = fft_field[:, :, : fft_field.shape[-1] // 2]
    freqs = np.fft.fftfreq(n_sample_points, d=dt)[: (int(n_sample_points / 2))]

    # We are only interested on the frequencies inside the source bandwidth
    mask = (freqs >= min_freq) * (freqs <= max_freq)
    freqs = freqs[mask]
    fft_field = fft_field[:, :, mask]

    return (freqs, fft_field)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot frequency domain field profile monitor."
    )
    parser.add_argument("--fname", help="Input h5 filename")
    parser.add_argument("--ename", help="Input h5 epsilon filename", default="")
    parser.add_argument(
        "--npoints",
        default=100,
        type=int,
        help="Number of frequency points to compute. Spacing is uniform.",
    )
    parser.add_argument(
        "--t0", default=0, type=float, help="Center of the apodization function."
    )
    parser.add_argument(
        "--tau",
        default=1e-10,
        type=float,
        help="Ramp up time of the apodization function.",
    )
    args = parser.parse_args()

    (field, metadata) = load_monitor(args.fname)

    apodize = lambda t: sigmoid_apodization(t, args.t0, args.tau)
    (freqs, fft_field) = transform_field(field, metadata, args.npoints, apodize)

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.1, bottom=0.25)

    cmap_fields = plt.cm.viridis

    # Check if epsilon argument supplied
    if args.ename != "":
        eps = load_struc(args.ename)

        # Choose colormap
        cmapbuf = plt.cm.Wistia

        # Get the colormap colors
        cmap_fields = cmapbuf(np.arange(cmapbuf.N))

        # Set alpha
        cmap_fields[:, -1] = np.linspace(0, 1, cmapbuf.N)

        # Create new colormap
        cmap_fields = ListedColormap(cmap_fields)

        # Step colormap for epsilon structure
        cmap_structure = ListedColormap(["navy", "black"])
        eps_plot = plt.imshow(
            eps,
            **default_eps_parameters,
            cmap=cmap_structure,
            extent=[
                metadata["xmin"],
                metadata["xmax"],
                metadata["ymin"],
                metadata["ymax"],
            ]
        )

    # Plot Results
    central_freq = (metadata["max_freq"] + metadata["min_freq"]) * 0.5
    field_plot = plt.imshow(
        fft_field[:, ::-1, find_nearest(freqs, central_freq)].T,
        cmap=cmap_fields,
        extent=[metadata["xmin"], metadata["xmax"], metadata["ymin"], metadata["ymax"]],
    )
    plt.colorbar()

    axcolor = "lightgoldenrodyellow"
    axfreq = plt.axes([0.1, 0.1, 0.65, 0.03], facecolor=axcolor)
    sfreq = Slider(
        axfreq,
        "Freq",
        np.min(freqs),
        np.max(freqs),
        valinit=central_freq,
        valfmt="%0.3f",
    )

    def update(val):
        freq = sfreq.val
        field_plot.set_data(fft_field[:, ::-1, find_nearest(freqs, freq)].T)
        field_plot.autoscale()
        fig.canvas.draw_idle()

    sfreq.on_changed(update)

    plt.show()
