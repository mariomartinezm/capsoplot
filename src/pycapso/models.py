import numpy as np


def get_mean_field(num_iter=100, psi0=1, phi0=0.01, alpha=0.1,
                   ey=1, ry=1, ez=1, rz=1,
                   a=-1, b=1, d=1, e=0):
    """Gets the mean field model approximation of the CAPSO model.

    Calculates the densities of preys and predators, using the mean field
    difference equations that describe the following interactions between
    individuals of both populations: intraspecific competition, reproduction of
    predators, death of predators, death of preys and reproduction of
    predators.

    Args:
        num_iter (int): the number of iterations to simulate.
        psi0 (float): the initial density of preys.
        phi0 (float): the initial density of predators.
        alpha (float) the intraspecific competition coefficient.
        ey (float): the reproductive capacity of preys.
        ry (int): the radius of the preys' reproduction neighborhood.
        ez (float): the reproductive capacity of predators.
        rz (int): the radius of the predators' reproduction neighborhood.
        a (float): the coefficient of the line that defines predator mortality.
        b (float): the intercept of the line that defines predator mortality.
        d (float): the coefficient of the line that defines prey mortality.
        e (float): the intercept of the line that defines prey mortality.

    Returns:
        A two dimensional numpy array, where the first column is a vector of
        prey densities and the second column is a vector of predators
        densities.
    """
    # Initialize data arrays
    index_set = np.arange(0, num_iter)
    psi = np.zeros(len(index_set))
    phi = np.zeros(len(index_set))

    # Calculate the probability of an 'event' occurring in the neighborhood of
    # a prey or predator
    card_mry = (2 * ry + 1) ** 2 - 1
    py = 1 / card_mry
    card_mrz = (2 * rz + 1) ** 2 - 1
    pz = 1 / card_mrz

    # Initialize densities and simulate the model
    psi[0] = psi0
    phi[0] = phi0

    for t in index_set[1:]:
        # Intraspecific competition
        psi_ic = psi[t - 1] - alpha * psi[t - 1] ** 2

        # Reproduction of predators
        num_preds = card_mrz * phi[t - 1]
        max_num_of_births = ez * num_preds
        phi_r = phi[t - 1] + (1 - phi[t - 1]) * (1 - (1 - pz) **
                                                 max_num_of_births)
        # Death of predators
        phi[t] = phi_r - (b + a * psi_ic) * phi_r

        # Death of preys
        psi_d = psi_ic - (e + d * phi[t]) * psi_ic

        # Reprodution of preys
        num_preys = card_mry * psi_d
        max_num_of_births = ey * num_preys
        psi[t] = psi_d + (1 - psi_d) * (1 - (1 - py) ** max_num_of_births)

    return np.array([psi, phi], ndmin=2).T
