import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import rfft, rfftfreq


def plot_prey_pred_data(data,
                        prey_label='Preys', pred_label='Predators',
                        prey_color='g', pred_color='r',
                        prey_linewidth=1.0, pred_linewidth=1.0,
                        prey_style='-', pred_style='-',
                        prey_marker='', pred_marker=''):
    """ A matplotlib.pyplot.plot() wrapper to quickly plot prey-predator data.

    Args:
        data (numpy array): the data to be plotted.
        prey_label (str): the label for the preys' data.
        pred_label (str): the label for the predators' data.
        prey_color (str): the color string for the preys' data.
        pred_color (str): the color string for the predators' data.
        prey_linewidth (float): the linewidth for the preys.
        pred_linewidth (float): the linewidth for the predators.
        prey_style (str): the string specifying the line style for the preys.
        pred_style (str): the string specifying the line style for the
            predators.
        prey_marker (str): the marker for the preys' data.
        pred_marker (str): the marker for the predators' data.
    """

    plt.plot(data[:, 0], antialiased=True, label=prey_label,
             color=prey_color, linewidth=prey_linewidth, linestyle=prey_style,
             marker=prey_marker)
    plt.plot(data[:, 1], antialiased=True, label=pred_label,
             color=pred_color, linewidth=pred_linewidth, linestyle=pred_style,
             marker=pred_marker)

    plt.xlabel('Time (Seasons)')
    plt.ylabel('Population density')

    plt.legend(loc='best')


def plot_fourier_spectrum(ts, sampling_freq=1):
    """ Plots the Fourier spectra of data in a Pandas series.

    Args:
        ts (Pandas.Series): The Pandas series that contains the data to be
        processed.
        sampling_freq (float): The sampling frequency of the data.
    """

    N = ts.size

    # Since the fft() function expects a numpy array we use the values property
    # of the series
    zero_mean_data = (ts - ts.mean()).values

    yf = rfft(zero_mean_data)
    xf = rfftfreq(N, 1 / sampling_freq)

    plt.plot(xf, np.abs(yf), 'k-', antialiased=True, linewidth=1.0)
    plt.xlabel("Frequency")
