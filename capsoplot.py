"""
    This module contains some plot utilites used to analyse the behavior of the
    Ca-Pso model.
"""

from __future__ import division
from pylab import figure, grid, rc, xlabel, ylabel, plot, gca, tick_params, legend, show
from scipy import loadtxt, arange, mean, fft, zeros, log
from scipy.fftpack import fftshift

AXIS_LABEL_SIZE = 25
TICKS_LABEL_SIZE = 20


def _set_font():
    rc('font', family='serif', serif='Helvetica')
    rc('text', usetex=True)


def _setup_grid_and_axes(label_x, label_y):
    grid(True)

    # Set axes labels
    xlabel(label_x, fontsize=AXIS_LABEL_SIZE)
    ylabel(label_y, fontsize=AXIS_LABEL_SIZE)

    # Set the axis ticks
    tick_params(axis='both', which='major', labelsize=TICKS_LABEL_SIZE)


def plot_time_series(file_name, tmin=-1, tmax=-1):
    """
    Plot the time series of a CaPso results file.

    Args:
        file_name (str)   -- the text file containing the results.

    Keyword arguments:
        tmin=-1 (int)     -- the minimal endpoint of the time interval to plot.
        tmax=-1 (int)     -- the maximum endpoint of the time interval to plot.

    """
    fig = figure(1)

    _set_font()

    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    ax1.set_xlabel('Time (Seasons)', fontsize=AXIS_LABEL_SIZE)
    ax1.set_ylabel('Preys', fontsize=AXIS_LABEL_SIZE)
    ax1.tick_params(axis='both', which='major', labelsize=TICKS_LABEL_SIZE)

    ax2 = ax1.twinx()
    ax2.set_ylabel('Predators', fontsize=AXIS_LABEL_SIZE)
    ax2.tick_params(axis='both', which='major', labelsize=TICKS_LABEL_SIZE)

    # load the results file
    index, preys, predators = loadtxt(file_name, unpack=True)

    # plot the prey's data
    if tmin != -1 and tmax != -1:
        preys_plot, = ax1.plot(index[tmin:tmax], preys[tmin:tmax], linewidth=1.5)
    else:
        preys_plot, = ax1.plot(index, preys, linewidth=1.5)

    preys_plot.set_antialiased(True)
    preys_plot.set_color('g')

    # plot the predator's data
    if tmin != -1 and tmax != -1:
        predators_plot, = ax2.plot(index[tmin:tmax], predators[tmin:tmax], linewidth=1.5)
    else:
        predators_plot, = ax2.plot(index, predators, linewidth=1.5)
    predators_plot.set_antialiased(True)
    predators_plot.set_color('r')

    legend([preys_plot, predators_plot], ['Preys', 'Predators'], loc=2, shadow=True)

    # Show the plot
    show()


def plot_time_series_normalized(file_name, tmin=-1, tmax=-1):
    """
    Plot the time series (normalized) of a CaPso results file.

    Args:
        file_name (str)   -- the text file containing the results.

    Keyword arguments:
        tmin=-1 (int)     -- the minimal endpoint of the time interval to plot.
        tmax=-1 (int)     -- the maximum endpoint of the time interval to plot.

    """
    size = 131072

    figure(1)

    _set_font()

    _setup_grid_and_axes('Time (Seasons)', 'Density')

    # load the results file
    index, preys, predators = loadtxt(file_name, unpack=True)

    # plot the prey's data
    if tmin != -1 and tmax != -1:
        plot(index[tmin:tmax], preys[tmin:tmax] / size, linewidth=1.5, label='Preys')
        plot(index[tmin:tmax], predators[tmin:tmax] / size, linewidth=1.5, label='Predators')
    else:
        plot(index, preys / size, linewidth=1.5, label='Preys')
        plot(index, predators / size, linewidth=1.5, label='Predators')

    legend()

    # Show the plot
    show()


def plot_fourier_spectra(file_name):
    # load data file
    index, preys, predators = loadtxt(file_name, unpack=True)

    N = preys.size

    f = arange(-N / 2, N / 2) / N

    zero_mean_data = preys - mean(preys)

    transform = fft(zero_mean_data)

    transform_scaled = transform / N

    F = abs(fftshift(transform_scaled))

    figure(1, (9, 8))

    _set_font()

    _setup_grid_and_axes(r'$\omega / 2 \pi$', '')

    # plot using a solid line
    plot(f, F, 'k-', antialiased=True, linewidth=1.0)
    x_axis = gca()
    x_axis.set_xlim([0, f.max()])

    # Show the plot
    show()


def plot_phase_plot(file_name, tmin=-1, tmax=-1):
    # Load data file
    index, preys, predators = loadtxt(file_name, unpack=True)

    figure(1, (9.0, 7))

    grid(True)

    _set_font()

    _setup_grid_and_axes('Preys', 'Predators')

    # Plot the data
    if tmin != -1 and tmax != -1:
        plot(preys[tmin:tmax] / 131072, predators[tmin:tmax] / 131072, 'k-', antialiased=True)
    else:
        plot(preys / 131072, predators / 131072, 'k-', antialiased=True)

    # Show the plot
    show()


def plot_mf_intraspecific(N=100, y0=1, alpha=0.5, z=0, data_file=''):
    # Initialize data arrays
    index_set = arange(0, N + 1)
    Y = zeros(len(index_set))
    preys = []

    if data_file != '':
        # Obtain data from data file
        index, preys, predators = loadtxt(data_file, unpack=True)
        preys = preys / 131072
        predators = predators / 131072

    # initialize densities
    Y[0] = y0

    # Calculate densities
    for t in index_set[1:]:
        Y[t] = Y[t - 1] - alpha * Y[t - 1] ** 2 + alpha * z * Y[t - 1]

    # Print raw data
    if data_file != '':
        print preys
    print Y

    # Setup the plot
    figure(1, (9.0, 7.0))

    _setup_grid_and_axes('t (seasons)', 'Population density')

    # Plot the data
    if data_file != '':
        plot(index_set, preys[0:N + 1], 'ro-', antialiased=True, label='Simulaton')
    plot(index_set, Y, 'bo-', antialiased=True, label='Mean field')

    legend()

    # Show the plot
    show()


def plot_mf_prey_reproduction(N=100, psi0=0.001, data_file=''):
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))

    preys = []
    predators = []

    if data_file != '':
        # Obtain data from file
        index, preys, predators = loadtxt(data_file, unpack=True)
        preys = preys / 131072
        predators = predators / 131072

    # Initialize densities
    psi[0] = psi0

    # Calculate densities
    for t in index_set[1:]:
        psi[t] = 2 * psi[t - 1] - psi[t - 1] ** 2

    print preys
    print psi

    # Setup the plot
    figure(1, (9.0, 7.0))

    _setup_grid_and_axes('t (seasons)', 'Population density')

    # Plot the data
    if data_file != '':
        plot(index_set, preys[0:N + 1], 'ro-', antialiased=True, label='Simulation')

    plot(index_set, psi, 'bo-', antialiased=True, label='Mean field')

    legend()

    # Show the plot
    show()


def plot_mf_only_preys(N=1000, tmin=-1, tmax=-1, psi0=0.0001, alpha=0.5):
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))

    # Initialize densities
    psi[0] = psi0

    # Calculate densities
    for t in index_set[1:]:
        psi[t] = psi[t - 1] - alpha * psi[t - 1] ** 2 + (1 - psi[t - 1]) * psi[t - 1]

    # Set up the plot
    figure(1)

    _setup_grid_and_axes('t (Seasons)', 'Population density')

    # Plot the data
    if tmin != -1 and tmax != -1:
        plot(index_set[tmin:tmax], psi[tmin:tmax], 'k.-', antialiased=True)
    else:
        plot(index_set, psi, 'k.-', antialiased=True)

    # Show the plot
    show()


def plot_mf_only_preys_coupled(N=1000, tmin=-1, tmax=-1, psi0=0.0001, alpha=0.5):
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))

    # Initialize densities
    psi[0] = psi0

    # Calculate densities
    for t in index_set[1:]:
        aux = psi[t - 1] - alpha * psi[t - 1] ** 2
        psi[t] = aux + (1 - aux) * aux

    # Set up the plot
    figure(1)

    _setup_grid_and_axes('t (Seasons)', 'Population density')

    # Plot the data
    if tmin != -1 and tmax != -1:
        plot(index_set[tmin:tmax], psi[tmin:tmax], 'k.-', antialiased=True)
    else:
        plot(index_set, psi, 'k.-', antialiased=True)

    # Show the plot
    show()


def plot_mf_prey_reproduction_2(N=100, psi0=0.001, ey=1, data_file=''):
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))

    preys = []
    predators = []

    if data_file != '':
        # Obtain data from file
        index, preys, predators = loadtxt(data_file, unpack=True)
        preys = preys / 131072
        predators = predators / 131072

    # Initialize densities
    psi[0] = psi0

    # Calculate densities
    for t in index_set[1:]:
        psi[t] = psi[t - 1] + (1 - psi[t - 1]) * ey * psi[t - 1]

    print preys
    print psi

    # Setup the plot
    figure(1, (9.0, 7.0))

    _setup_grid_and_axes('t (seasons)', 'Population density')

    # Plot the data
    if data_file != '':
        plot(index_set, preys[0:N + 1], 'ro-', antialiased=True, label='Simulation')

    plot(index_set, psi, 'bo-', antialiased=True, label='Mean field')

    legend()

    # Show the plot
    show()


def plot_mf_minimal(N=100, psi0=1, phi0=0.01, data_file=''):
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))
    phi = zeros(len(index_set))

    preys = []
    predators = []

    if data_file != '':
        # Obtain data from file
        index, preys, predators = loadtxt(data_file, unpack=True)
        preys = preys / 131072
        predators = predators / 131072

    # Initialize densities
    psi[0] = psi0
    phi[0] = phi0

    # Calculate densities
    for t in index_set[1:]:
        psi[t] = psi[t - 1] - phi[t - 1] * psi[t - 1]
        phi[t] = phi[t - 1] + (1 - phi[t - 1]) * phi[t - 1] - (1 - psi[t - 1]) * phi[t - 1]

    # Setup the plot
    figure(1)

    _setup_grid_and_axes('t (seasons)', 'Population density')

    # Plot the data
    if(data_file != ''):
        plot(index_set, preys[0:N + 1], 'co-', antialiased=True, label='Sim preys')
        plot(index_set, predators[0:N + 1], 'mo-', antialiased=True, label='Sim predators')

    plot(index_set, psi, 'go-', antialiased=True, label='Mf preys')
    plot(index_set, phi, 'ro-', antialiased=True, label='Mf predators')

    legend()

    # Show the plot
    show()


def plot_mf(N=100, psi0=1, phi0=0.01, alpha=0.1, ey=1, ez=1, data_file=''):
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))
    phi = zeros(len(index_set))

    preys = []
    predators = []

    if data_file != '':
        # Obtain data from file
        index, preys, predators = loadtxt(data_file, unpack=True)
        preys = preys / 131072
        predators = predators / 131072

    # Initialize densities
    psi[0] = psi0
    phi[0] = phi0

    # Calculate densities
    for t in index_set[1:]:
        psi[t] = psi[t - 1] + (1 - psi[t - 1]) * ey * psi[t - 1] - phi[t - 1] * psi[t - 1] - alpha * psi[t - 1] ** 2
        phi[t] = phi[t - 1] + (1 - phi[t - 1]) * ez * phi[t - 1] - (1 - psi[t - 1]) * phi[t - 1]

    # Setup the plot
    figure(1)

    _setup_grid_and_axes('t (seasons)', 'Population density')

    # Plot the data
    if(data_file != ''):
        plot(index_set, preys[0:N + 1], 'c-', antialiased=True, label='Sim preys')
        plot(index_set, predators[0:N + 1], 'm-', antialiased=True, label='Sim predators')

    plot(index_set, psi, 'g-', antialiased=True, label='Mf preys')
    plot(index_set, phi, 'r-', antialiased=True, label='Mf predators')

    legend()

    # Show the plot
    show()


def plot_mf_coupled(N=100, psi0=1, phi0=0.01, alpha=0.1, ey=1, ez=1, data_file=''):
# Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))
    phi = zeros(len(index_set))

    preys = []
    predators = []

    if data_file != '':
        # Obtain data from file
        index, preys, predators = loadtxt(data_file, unpack=True)
        preys = preys / 131072
        predators = predators / 131072

    # Initialize densities
    psi[0] = psi0
    phi[0] = phi0

    # Calculate densities
    for t in index_set[1:]:
        psi_ic = psi[t - 1] - alpha * psi[t - 1] ** 2
        phi_rz = phi[t - 1] + (1 - phi[t - 1]) * ez * phi[t - 1]
        phi[t] = phi_rz - (1 - psi_ic) * phi_rz
        psi_dy = psi_ic - phi[t] * psi_ic
        psi[t] = psi_dy + (1 - psi_dy) * ey * psi_dy

    # Setup the plot
    figure(1)

    _setup_grid_and_axes('t (seasons)', 'Population density')

    # Plot the data
    if(data_file != ''):
        plot(index_set, preys[0:N + 1], 'c-', antialiased=True, label='Sim preys')
        plot(index_set, predators[0:N + 1], 'm-', antialiased=True, label='Sim predators')

    plot(index_set, psi, 'g-', antialiased=True, label='Mf preys')
    plot(index_set, phi, 'r-', antialiased=True, label='Mf predators')

    legend()

    # Show the plot
    show()


def plot_mf_coupled_phase(N=100, psi0=1, phi0=0.01, alpha=0.1, ey=1, ez=1):
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))
    phi = zeros(len(index_set))

    # Initialize densities
    psi[0] = psi0
    phi[0] = phi0

    # Calculate densities
    for t in index_set[1:]:
        psi_ic = psi[t - 1] - alpha * psi[t - 1] ** 2
        phi_rz = phi[t - 1] + (1 - phi[t - 1]) * ez * phi[t - 1]
        phi[t] = phi_rz - (1 - psi_ic) * phi_rz
        psi_dy = psi_ic - phi[t] * psi_ic
        psi[t] = psi_dy + (1 - psi_dy) * ey * psi_dy

    # Setup the plot
    figure(1)

    _setup_grid_and_axes('t (Seasons)', 'Population density')

    # Plot the data
    plot(psi, phi, 'k-', antialiased=True)

    # Show the plot
    show()


def plot_mf_phase(N=100, tmin=-1, tmax=-1, psi0=1, phi0=0.01, alpha=0.1, ey=1, ez=1):
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))
    phi = zeros(len(index_set))

    # Initialize densities
    psi[0] = psi0
    phi[0] = phi0

    # Calculate densities
    for t in index_set[1:]:
        psi[t] = psi[t - 1] + (1 - psi[t - 1]) * ey * psi[t - 1] - phi[t - 1] * psi[t - 1] - alpha * psi[t - 1] ** 2
        phi[t] = phi[t - 1] + (1 - phi[t - 1]) * ez * phi[t - 1] - (1 - psi[t - 1]) * phi[t - 1] - phi[t - 1]

    # Setup the plot
    figure(1)

    _setup_grid_and_axes('Preys', 'Predators')

    # Plot the data
    if tmin != -1 and tmax != -1:
        plot(psi[tmin:tmax], phi[tmin:tmax], 'k-', antialiased=True)
    else:
        plot(psi, phi, 'k-', antialiased=True)

    # Show the plot
    show()


def plot_mf_seck(N=1000, tmin=-1, tmax=-1, psi0=0.0001, alpha=0.5, rc=1, rrp=1):
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))

    # Initialize densities
    psi[0] = psi0

    # Initialize extra parameters
    cmc = 1 / ((2 * rc + 1) ** 2)
    aux = 0.6
    w = 1 - (aux ** rrp)

    # Calculate densities
    for t in index_set[1:]:
        dv = (1 / (1 - cmc)) * psi[t - 1] + 1 - (1 / (1 - cmc))
        #psi[t] = psi[t - 1] - alpha * psi[t - 1] * dv
        #psi[t] = psi[t - 1] + (w / log(t + 1)) * (1 - psi[t - 1]) * psi[t - 1]
        psi[t] = psi[t - 1] - alpha * psi[t - 1] * dv + (w / log(t + 1)) * (1 - psi[t - 1]) * psi[t - 1]

    # Set up the plot
    figure(1)

    _setup_grid_and_axes('t (Seasons)', 'Population density')

    # Plot the data
    if tmin != -1 and tmax != -1:
        plot(index_set[tmin:tmax], psi[tmin:tmax], 'k.-', antialiased=True)
    else:
        plot(index_set, psi, 'k.-', antialiased=True)

    # Show the plot
    show()


def plot_mf_seck_coupled(N=1000, tmin=-1, tmax=-1, psi0=0.0001, alpha=0.5, rc=1, rrp=1):
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))

    # Initialize densities
    psi[0] = psi0

    # Initialize extra parameters
    cmc = 1 / ((2 * rc + 1) ** 2)
    aux = 0.6
    w = 1 - (aux ** rrp)

    # Calculate densities
    for t in index_set[1:]:
        dv = (1 / (1 - cmc)) * psi[t - 1] + 1 - (1 / (1 - cmc))
        temp = psi[t - 1] - alpha * psi[t - 1] * dv
        psi[t] = temp + (w / log(t + 1)) * (1 - temp) * temp

    # Set up the plot
    figure(1)

    _setup_grid_and_axes('t (Seasons)', 'Population density')

    # Plot the data
    if tmin != -1 and tmax != -1:
        plot(index_set[tmin:tmax], psi[tmin:tmax], 'k.-', antialiased=True)
    else:
        plot(index_set, psi, 'k.-', antialiased=True)

    # Show the plot
    show()
