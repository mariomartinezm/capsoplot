"""
    This module contains some plot utilites used to analyse the behavior of the
    Ca-Pso model.
"""

from __future__ import division
from os import listdir
from os.path import isfile, join
from pylab import figure, grid, rc, xlabel, ylabel, plot, gca, tick_params
from pylab import legend, show
from scipy import loadtxt, arange, mean, fft, zeros, polyfit, polyval, sqrt
from scipy.fftpack import fftshift
from scipy.optimize import leastsq

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


def get_average_column(path, column=0):
    """
    Get the index-based average column for a series of results files.

    Args:
        path(str): the path containing the results files.

    Kwargs:
        column (int): the column index in a results file.

    Returns:
        A numpy.ndarray containing the average values for the specified
        column-index of a series of results files.

    """
    files = [f for f in listdir(path) if isfile(join(path, f))
             and f.endswith(".txt")]

    col_seq = column,

    sum_col = loadtxt(join(path, files[0]), usecols=col_seq, unpack=True)

    for file in files[1:]:
        sum_col = sum_col + loadtxt(join(path, file), usecols=col_seq,
                                    unpack=True)

    return sum_col / len(files)


def plot_time_series(file_name, tmin=-1, tmax=-1):
    """
    Plot the time series of a CaPso results file.

    Args:
        file_name (str): the text file containing the results.

    Kwargs:
        tmin (int): the minimal endpoint of the time interval to plot.
        tmax (int): the maximum endpoint of the time interval to plot.

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
    index, preys, predators = loadtxt(file_name, usecols=(0, 1, 2),
                                      unpack=True)

    # plot the prey's data
    if tmin != -1 and tmax != -1:
        preys_plot, = ax1.plot(index[tmin:tmax], preys[tmin:tmax],
                               linewidth=1.5)
    else:
        preys_plot, = ax1.plot(index, preys, linewidth=1.5)

    preys_plot.set_antialiased(True)
    preys_plot.set_color('g')

    # plot the predator's data
    if tmin != -1 and tmax != -1:
        predators_plot, = ax2.plot(index[tmin:tmax], predators[tmin:tmax],
                                   linewidth=1.5)
    else:
        predators_plot, = ax2.plot(index, predators, linewidth=1.5)
    predators_plot.set_antialiased(True)
    predators_plot.set_color('r')

    legend([preys_plot, predators_plot], ['Preys', 'Predators'], loc=2,
           shadow=True)

    # Show the plot
    show()


def plot_time_series_normalized(file_name, width=512, height=256, tmin=-1,
                                tmax=-1, prey_label='Preys',
                                pred_label='Predators', prey_color='g',
                                pred_color='r', prey_style='-',
                                pred_style='-', prey_marker='',
                                pred_marker=''):
    """
    Plot the time series (normalized) of a CaPso results file.

    Args:
        file_name (str): the text file containing the data to plot.

    Kwargs:
        width (int): the width of the lattice used in the simulation.
        height (int): the height of the lattice used in the simulation.
        tmin (int): the minimal endpoint of the time interval to plot.
        tmax (int): the maximum endpoint of the time interval to plot.
        prey_label (str): the label for the preys' data.
        pred_label (str): the label for the predators' data.
        prey_color (str): the color string for the preys' data.
        pred_color (str): the color string for the predators' data.
        prey_style (str): the string specifying the line style for the preys.
        pred_style (str): the string specifying the line style for the
            predators.
        prey_marker (str): the marker for the preys' data.
        pred_marker (str): the marker for the predators' data.

    """
    size = width * height

    figure(1, (9, 8))

    _set_font()

    _setup_grid_and_axes('Time (Seasons)', 'Density')

    # load the results file
    index, preys, predators = loadtxt(file_name, usecols=(0, 1, 2),
                                      unpack=True)

    # plot the prey's data
    if tmin != -1 and tmax != -1:
        plot(index[tmin:tmax], preys[tmin:tmax] / size, label=prey_label,
             color=prey_color, linestyle=prey_style, marker=prey_marker)
        plot(index[tmin:tmax], predators[tmin:tmax] / size, label=pred_label,
             color=pred_color, linestyle=pred_style, marker=pred_marker)
    else:
        plot(index, preys / size, label=prey_label, color=prey_color,
             linestyle=prey_style, marker=prey_marker)
        plot(index, predators / size, label=pred_label, color=pred_color,
             linestyle=pred_style, marker=pred_marker)

    legend()

    # Show the plot
    show()


def plot_fourier_spectra(file_name):
    """
    Plot the the Fourier spectra of a CaPso results file.

    Args:
        file_name (str): the text file containing the results.

    """
    # load data file
    index, preys = loadtxt(file_name, usecols=(0, 1), unpack=True)

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


def plot_phase_plot(file_name, width=512, height=256, tmin=-1, tmax=-1,
                    label='Sim', color='k', style='-', marker=''):
    """
    Plot the the phase plot of a CaPso results file.

    Args:
        file_name (str): the text file containing the results.

    Kwargs:
        width (int): the width of the lattice used in the simulation.
        height (int): the height of the lattice used in the simulation.
        tmin (int): the minimal endpoint of the time interval to plot.
        tmax (int): the maximum endpoint of the time interval to plot.
        label (string): the label for the data.
        color (string): the string specifying the line color.
        style (string): the string specifying the line style.
        marker (string): the string specifying the marker for the data points.

    """
    size = width * height

    # Load data file
    index, preys, predators = loadtxt(file_name, usecols=(0, 1, 2),
                                      unpack=True)

    figure(1, (9.0, 8))

    grid(True)

    _set_font()

    _setup_grid_and_axes('Preys', 'Predators')

    # Plot the data
    if tmin != -1 and tmax != -1:
        plot(preys[tmin:tmax] / size, predators[tmin:tmax] / size,
             label=label, color=color, linestyle=style, marker=marker,
             antialiased=True)
    else:
        plot(preys / size, predators / size, label=label, color=color,
             linestyle=style, marker=marker, antialiased=True)

    legend()

    # Show the plot
    show()


def plot_birth_rate(file_name, use_prey_data=True, width=512, height=256,
                    label='Prey birth rate', color='k', style='.', marker='.'):
    """
    Plot the birth rate using a CaPso results file.

    Args:
        file_name (str): the text file containing the results.

    Kwargs:
        use_prey_data (bool): if true, the birth rate of preys is plotted,
        otherwise, the birth rate of predators is plotted.
        width (int): the width of the lattice used in the simulation.
        height (int): the height of the lattice used in the simulation.
        label (string): the label for the data.
        color (string): the string specifying the line color.
        style (string): the string specifying the line style.
        marker (string): the string specifying the marker for the data points.

    """
    if use_prey_data is True:
        pop_size, birth_rate = loadtxt(file_name, usecols=(3, 4), unpack=True)
    else:
        pop_size, birth_rate = loadtxt(file_name, usecols=(5, 6), unpack=True)

    size = width * height

    figure(1, (9.0, 7))

    grid(True)

    _set_font()

    _setup_grid_and_axes('Population density', 'Birth rate')

    # Plot the data
    plot(pop_size / size, birth_rate, label=label, color=color,
         linestyle=style, marker=marker, antialiased=True)

    legend()

    # Show the plot
    show()


def plot_leastsq_for_reproduction(pop, birth_rate, e, r,
                                  x_label='Density',
                                  y_label='Growth rate',
                                  fit_label='Fitted curve',
                                  birth_rate_label='Birth rate data',
                                  fit_color='b',
                                  birth_rate_color='g',
                                  fit_style='-',
                                  birth_rate_style='.',
                                  fit_marker='',
                                  birth_rate_marker='.'):
    p0 = [e, r]

    plsq = leastsq(_residuals, p0, args=(birth_rate, pop), full_output=1)

    print 'Solution found: {0}'.format(plsq[0])
    print 'Output flag: {0}'.format(plsq[4])
    print 'Output msg: {0}'.format(plsq[3])

    figure(1, (9, 8))

    _set_font()

    _setup_grid_and_axes(x_label, y_label)

    plot(pop, _peval(pop, plsq[0]), label=fit_label, color=fit_color,
         linestyle=fit_style, marker=fit_marker)
    plot(pop, birth_rate, '.', label=birth_rate_label, color=birth_rate_color,
         linestyle=birth_rate_style, marker=birth_rate_marker)

    legend()

    show()


def _residuals(params, y, x):
    e, r = params
    card = (2 * r + 1) ** 2 - 1
    p = 1 / card
    err = y - (1 - x) * (1 - (1 - p) ** (card * e * x))

    return err


def _peval(x, params):
    e = params[0]
    r = params[1]
    card = (2 * r + 1) ** 2 - 1
    p = 1 / card

    return (1 - x) * (1 - (1 - p) ** (card * e * x))


def plot_reg_for_predator_death(prey_data, pred_dp):
    n = len(prey_data)

    # polynomial regression
    (ar, br) = polyfit(prey_data, pred_dp, 1)
    reg = polyval([ar, br], prey_data)

    # compute the mean square error
    err = sqrt(sum((reg - pred_dp) ** 2) / n)

    print('Parameters: a = {0}, b = {1}'.format(ar, br))
    print('Mean square error = {0}'.format(err))

    # plots
    plot(prey_data, pred_dp, 'b.')
    plot(prey_data, reg, 'g-')
    plot(prey_data, 1 - prey_data, 'r-')
    legend(['Original data', 'Regression', 'Mean field'])


def plot_reg_for_prey_death(pred_data, prey_dp):
    n = len(pred_data)

    # polynomial regression
    cr = polyfit(pred_data, prey_dp, 1)
    reg = polyval(cr, pred_data)

    # compute the mean square error
    err = sqrt(sum((reg - prey_dp) ** 2) / n)

    print('Parameters: c = {0}'.format(cr))
    print('Mean square error = {0}'.format(err))

    # plots
    plot(pred_data, prey_dp, 'b.')
    plot(pred_data, reg, 'g-')
    plot(pred_data, pred_data, 'r-')
    legend(['Original data', 'Regression', 'Mean field'])


def plot_mf_intraspecific(N=100, y0=1, alpha=0.5, z=0, data_file=''):
    """
    Plot the mean field term for the intraspecific competition.

    Kwargs:
        N (int): the number of iterations to calculate.
        y0 (float): the initial density of preys.
        alpha (float): the intraspecific competition coefficient.
        z (float): an adjustement factor.
        data_file (str): an optional CaPso results file.

    """
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
        plot(index_set, preys[0:N + 1], 'ro-', antialiased=True,
             label='Simulaton')

    plot(index_set, Y, 'bo-', antialiased=True, label='Mean field')

    legend()

    # Show the plot
    show()


def plot_mf_prey_reproduction(N=100, psi0=0.001, ry=1, ey=1, data_file=''):
    """
    Plot the mean field term for the reproduction stage.

    Kwargs:
        N (int): the number of iterations to calculate.
        psi0 (float): the initial density of preys.
        ry (int): the radious of the reproduction neighborhood.
        ey (int): the reproductive capacity of the preys.
        data_file (str): an optional CaPso results file.

    """
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

    card_mry = (2 * ry + 1) ** 2 - 1
    py = 1 / card_mry

    # Calculate densities
    for t in index_set[1:]:
        number_of_preys = card_mry * psi[t - 1]
        number_of_events = ey * number_of_preys

        psi[t] = psi[t - 1] + (1 - psi[t - 1]) * (1 - (1 - py) **
                                                  number_of_events)

    # Setup the plot
    figure(1)

    _setup_grid_and_axes('t (seasons)', 'Population density')

    # Plot the data
    if data_file != '':
        plot(index_set, preys[0:N + 1], 'b-', antialiased=True,
             label='Simulation')

    plot(index_set, psi, 'g-', antialiased=True, label='Mf preys')

    legend()

    # Show the plot
    show()


def plot_mf_minimal(N=100, psi0=1, phi0=0.01, ez=1, rz=1,
                    prey_label='Mf preys', pred_label='Mf predators',
                    prey_color='g', pred_color='r', prey_style='-',
                    pred_style='-', prey_marker='', pred_marker=''):
    """
    Plot a mean field model that does not include intraspecific competiton nor
    prey reproduction.

    Kwargs:
        N (int): the number of iterations to calculate.
        psi0 (float): the initial density of preys.
        phi0 (float): the initial density of predators.
        ez (int): the reproductive capacity of predators.
        rz (int): the radius of the reproduction neighborhood.
        prey_label (str): the label for the preys' data.
        pred_label (str): the label for the predators' data.
        prey_color (str): the color string for the preys' data.
        pred_color (str): the color string for the predators' data.
        prey_style (str): the string specifying the line style for the preys.
        pred_style (str): the string specifying the line style for the
            predators.
        prey_marker (str): the marker for the preys' data.
        pred_marker (str): the marker for the predators' data.

    """
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))
    phi = zeros(len(index_set))

    # Initialize densities
    psi[0] = psi0
    phi[0] = phi0

    card_mrz = (2 * rz + 1) ** 2 - 1
    pz = 1 / card_mrz

    # Calculate densities
    for t in index_set[1:]:
        nprz = (1 - pz) ** (ez * card_mrz * phi[t - 1])
        phi[t] = psi[t - 1] * (1 + nprz * (phi[t - 1] - 1) - phi[t - 1])
        psi[t] = psi[t - 1] - psi[t - 1] * ((1 + nprz * (phi[t - 1] - 1)) -
                                            phi[t - 1])

    # Setup the plot
    figure(1)

    _set_font()

    _setup_grid_and_axes('t (seasons)', 'Population density')

    plot(index_set, psi, antialiased=True, label=prey_label, color=prey_color,
         linestyle=prey_style, marker=prey_marker)
    plot(index_set, phi, antialiased=True, label=pred_label, color=pred_color,
         linestyle=pred_style, marker=pred_marker)

    legend()

    # Show the plot
    show()


def plot_mf_minimal_coupled(N=100, psi0=1, phi0=0.01, ez=1, rz=1,
                            prey_label='Mf preys', pred_label='Mf predators',
                            prey_color='g', pred_color='r', prey_style='-',
                            pred_style='-', prey_marker='', pred_marker=''):
    """
    Plot a mean field model that does not include intraspecific competiton nor
    prey reproduction (coupled version).

    Kwargs:
        N (int): the number of iterations to calculate.
        psi0 (float): the initial density of preys.
        phi0 (float): the initial density of predators.
        ez (int): the reproductive capacity of predators.
        rz (int): the radius of the reproduction neighborhood.
        prey_label (str): the label for the preys' data.
        pred_label (str): the label for the predators' data.
        prey_color (str): the color string for the preys' data.
        pred_color (str): the color string for the predators' data.
        prey_style (str): the string specifying the line style for the preys.
        pred_style (str): the string specifying the line style for the
            predators.
        prey_marker (str): the marker for the preys' data.
        pred_marker (str): the marker for the predators' data.

    """
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))
    phi = zeros(len(index_set))

    # Initialize densities
    psi[0] = psi0
    phi[0] = phi0

    card_mrz = (2 * rz + 1) ** 2 - 1
    pz = 1 / card_mrz

    # Calculate densities
    for t in index_set[1:]:
        # Reproduction of predators
        number_of_predators = card_mrz * phi[t - 1]
        number_of_events_predators = ez * number_of_predators
        phi_rz = phi[t - 1] + \
            (1 - phi[t - 1]) * (1 - (1 - pz) ** number_of_events_predators)
        # Death of predators
        phi[t] = phi_rz - phi[t - 1] - (1 - psi[t - 1]) * (phi_rz - phi[t - 1])
        # Death of preys
        psi[t] = psi[t - 1] - phi[t]

    # Setup the plot
    figure(1)

    _set_font()

    _setup_grid_and_axes('t (seasons)', 'Population density')

    plot(index_set, psi, antialiased=True, label=prey_label, color=prey_color,
         linestyle=prey_style, marker=prey_marker)
    plot(index_set, phi, antialiased=True, label=pred_label, color=pred_color,
         linestyle=pred_style, marker=pred_marker)

    legend()

    # Show the plot
    show()


def plot_mf(N=100, psi0=1, phi0=0.01, alpha=0.1, ey=1, ry=1, ez=1, rz=1,
            data_file=''):
    """
    Plot the mean field model approximation.

    Kwargs:
        N (int): the number of iterations to calculate.
        psi0 (float): the initial density of preys.
        phi0 (float): the initial density of predators.
        alpha (float): the intraspecific competition coefficient.
        ey (int): the reproductive capacity of preys.
        ry (int): the radius of the preys' reproduction neighborhood.
        ez (int): the reproductive capacity of predators.
        rz (int): the radius of the predators' reproduction neighborhood.
        data_file (str): an optional CaPso results file for comparison.

    """
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

    card_mry = (2 * ry + 1) ** 2 - 1
    py = 1 / card_mry

    card_mrz = (2 * rz + 1) ** 2 - 1
    pz = 1 / card_mrz

    # Calculate densities
    for t in index_set[1:]:
        ic = psi[t - 1] - alpha * psi[t - 1] ** 2
        pnrz = (1 - pz) ** (ez * phi[t - 1] * card_mrz)

        phi[t] = (1 + pnrz * (phi[t - 1] - 1)) * ic

        psi[t] = 1 + ((1 - py) ** (-ey * card_mry * ic * (phi[t - 1] - 1) *
                                   pnrz)) * (-ic * (phi[t - 1] - 1) * pnrz - 1)

    # Setup the plot
    figure(1)

    _set_font()

    _setup_grid_and_axes('t (seasons)', 'Population density')

    # Plot the data
    if(data_file != ''):
        plot(index_set, preys[0:N + 1], 'c-', antialiased=True,
             label='Sim preys')
        plot(index_set, predators[0:N + 1], 'm-', antialiased=True,
             label='Sim predators')

    plot(index_set, psi, 'g-', antialiased=True, label='Mf preys')
    plot(index_set, phi, 'r-', antialiased=True, label='Mf predators')

    legend()

    # Show the plot
    show()


def plot_mf_coupled(N=100, psi0=1, phi0=0.01, alpha=0.1, ey=1, ry=1, ez=1,
                    rz=1, a=-1, b=1, c=1, d=0, prey_label='Mf preys',
                    pred_label='Mf predators', prey_color='g', pred_color='r',
                    prey_style='-', pred_style='-', prey_marker='',
                    pred_marker=''):
    """
    Plot the mean field model approximation, coupled version.

    Kwargs:
        N (int): the number of iterations to calculate.
        psi0 (float): the initial density of preys.
        phi0 (float): the initial density of predators.
        alpha (float) the intraspecific competition coefficient.
        ey (float): the reproductive capacity of preys.
        ry (int): the radius of the preys' reproduction neighborhood.
        ez (float): the reproductive capacity of predators.
        rz (int): the radius of the predators' reproduction neighborhood.
        a (float): the ordinate of the line that defines predator mortality.
        b (float): the slope of the line that defines predator mortality.
        c (float): the slope of the line that defines prey mortality.
        prey_label (str): the label for the preys' data.
        pred_label (str): the label for the predators' data.
        prey_color (str): the color string for the preys' data.
        pred_color (str): the color string for the predators' data.
        prey_style (str): the string specifying the line style for the preys.
        pred_style (str): the string specifying the line style for the
            predators.
        prey_marker (str): the marker for the preys' data.
        pred_marker (str): the marker for the predators' data.

    """
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))
    phi = zeros(len(index_set))

    # Initialize densities
    psi[0] = psi0
    phi[0] = phi0

    card_mry = (2 * ry + 1) ** 2 - 1
    py = 1 / card_mry

    card_mrz = (2 * rz + 1) ** 2 - 1
    pz = 1 / card_mrz

    # Calculate densities
    for t in index_set[1:]:
        # Intraspecific competition
        psi_ic = psi[t - 1] - alpha * psi[t - 1] ** 2
        # Reproduction of predators
        number_of_predators = card_mrz * phi[t - 1]
        number_of_events_predators = ez * number_of_predators
        phi_rz = phi[t - 1] + \
            (1 - phi[t - 1]) * (1 - (1 - pz) ** number_of_events_predators)
        # Death of predators
        phi[t] = phi_rz - (b + a * psi_ic) * phi_rz
        # Death of preys
        psi_dy = psi_ic - (d + c * phi[t]) * psi_ic
        # Reprodution of preys
        number_of_preys = card_mry * psi_dy
        number_of_events = ey * number_of_preys
        psi[t] = psi_dy + (1 - psi_dy) * (1 - (1 - py) ** number_of_events)

    # Setup the plot
    figure(1, (9, 8))

    _set_font()

    _setup_grid_and_axes('t (seasons)', 'Population density')

    plot(index_set, psi, antialiased=True, label=prey_label, color=prey_color,
         linestyle=prey_style, marker=prey_marker)
    plot(index_set, phi, antialiased=True, label=pred_label, color=pred_color,
         linestyle=pred_style, marker=pred_marker)

    legend()

    # Show the plot
    show()


def plot_mf_coupled_phase(N=100, psi0=1, phi0=0.01, alpha=0.1, ey=1, ry=1,
                          ez=1, rz=1, label='Mf', color='k', style='-',
                          marker=''):
    """
    Plot the phase plot of the mean field model approximation, coupled version.

    Kwargs:
        N (int): the number of iterations to calculate.
        psi0 (float): the initial density of preys.
        phi0 (float): the initial density of predators.
        alpha (float): the intraspecific competition coefficient.
        ey (float): the reproductive capacity of preys.
        ry (int): the radius of the preys' reproduction neighborhood.
        ez (float): the reproductive capacity of predators.
        rz (int): the radius of the predators' reproduction neighborhood.
        label (string): the label for the data.
        color (string): the string specifying the line color.
        style (string): the string specifying the line style.
        marker (string): the string specifying the marker for the data points.

    """
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))
    phi = zeros(len(index_set))

    # Initialize densities
    psi[0] = psi0
    phi[0] = phi0

    card_mry = (2 * ry + 1) ** 2 - 1
    py = 1 / card_mry

    card_mrz = (2 * rz + 1) ** 2 - 1
    pz = 1 / card_mrz

    # Calculate densities
    for t in index_set[1:]:
        # Intraspecific competition
        psi_ic = psi[t - 1] - alpha * psi[t - 1] ** 2
        # Reproduction of predators
        number_of_predators = card_mrz * phi[t - 1]
        number_of_events_predators = ez * number_of_predators
        phi_rz = phi[t - 1] + \
            (1 - phi[t - 1]) * (1 - (1 - pz) ** number_of_events_predators)
        # Death of predators
        phi[t] = phi_rz - (1 - psi_ic) * phi_rz
        # Death of preys
        psi_dy = psi_ic - phi[t]
        # Reprodution of preys
        number_of_preys = card_mry * psi_dy
        number_of_events = ey * number_of_preys
        psi[t] = psi_dy + (1 - psi_dy) * (1 - (1 - py) ** number_of_events)

    # Setup the plot
    figure(1, (9, 8))

    _set_font()

    _setup_grid_and_axes('Preys', 'Predators')

    # Plot the data
    plot(psi, phi, 'k-', label=label, color=color, linestyle=style,
         marker=marker, antialiased=True)

    legend()

    # Show the plot
    show()


def plot_mf_phase(N=100, tmin=-1, tmax=-1, psi0=1, phi0=0.01, alpha=0.1,
                  ey=1, ez=1):
    """
    Plot the phase plot of the mean field model approximation.

    Kwargs:
        N (int): the number of iterations to calculate.
        tmin (int): the minimal endpoint of the time interval to plot.
        tmax (int): the maximum endpoint of the time interval to plot.
        psi0  (float): the initial density of preys.
        phi0 (float): the initial density of predators.
        alph (float): the intraspecific competition coefficient.
        ey (float): the reproductive capacity of preys.
        ez (float): the reproductive capacity of predators.

    """
    # Initialize data arrays
    index_set = arange(0, N + 1)
    psi = zeros(len(index_set))
    phi = zeros(len(index_set))

    # Initialize densities
    psi[0] = psi0
    phi[0] = phi0

    # Calculate densities
    for t in index_set[1:]:
        psi[t] = psi[t - 1] + (1 - psi[t - 1]) * ey * psi[t - 1] - \
            phi[t - 1] * psi[t - 1] - alpha * psi[t - 1] ** 2
        phi[t] = phi[t - 1] + (1 - phi[t - 1]) * ez * phi[t - 1] - \
            (1 - psi[t - 1]) * phi[t - 1] - phi[t - 1]

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
