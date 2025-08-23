#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Utility for evaluation vs. measurements of Irvine et al (1997)."""

from __future__ import division
from __future__ import unicode_literals

import os
import matplotlib 
matplotlib.rcParams["backend"] = "Agg"
from matplotlib import pyplot as plt
import numpy as np
from math import log, sqrt

REF_PROFILE = 'mast1'
PROFILE_NAMES = ('mast1', 'mast2', 'mast3', 'mast4')
TREE_HEIGHT = 7.5
REF_HEIGHT = 2 * TREE_HEIGHT
Z0 = 0.06
KARMAN = 0.41


def calculate_TI(sigmau, sigmav, sigmaw, u):
    """Calculate turbulence intensity from measurements of speed std.dev."""
    return np.sqrt(
        np.power(sigmau, 2) + np.power(sigmav, 2) + np.power(sigmaw, 2)
    ) / u


def ustar(uref, zref, z0):
    """Calculate friction velocity for a neutral ABL logartithmic profile."""
    return KARMAN * uref / log(zref / z0)


def read_U(filename):
    """ read z and U from vertical profile."""
    U = np.genfromtxt(filename, dtype=float)
    return (U[:, 0], U[:, 1:])


def read_epsilon_k_p(filename):
    """read epsilon, k p from vertical profile."""

    epsilon_k_p = np.genfromtxt(filename, dtype=float)
    z = epsilon_k_p[:, 0]
    epsilon = epsilon_k_p[:, 1]
    k = epsilon_k_p[:, 2]
    p = epsilon_k_p[:, 3]
    return z, epsilon, k, p


def main():

    times = map(int, os.listdir('postProcessing/verticalProfiles'))
    latestTime = max(times)
    fig1, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(
        4, 2, sharey='row', sharex='col'
    )

    filenameV = 'postProcessing/verticalProfiles/{time}/{profile}_vector.xy'.format(
        time=latestTime,
        profile=REF_PROFILE
    )
    filenameS = 'postProcessing/verticalProfiles/{time}/{profile}_scalar.xy'.format(
        time=latestTime,
        profile=REF_PROFILE
    )

    measurements = np.genfromtxt(
        'constant/experimentalData/measurements_irvine_1997.csv', delimiter=';',
        names=True, skip_header=1
    )

    z, Uref = read_U(filenameV)

    # model reference conditions (at mast 1, 2h)
    Uref_mag_model = np.linalg.norm(Uref, axis=1)
    ref_height_index = (np.abs(z - REF_HEIGHT)).argmin()
    uref_model = Uref_mag_model[ref_height_index]

    # measurements reference conditions (at mast 1, 2h)
    uref_measurement_norm = np.extract(
        (measurements['Mast'] == 1) & (measurements['z_n'] == 2),
        measurements
    )['U_n']

    ustar_ref_model = ustar(uref_model, REF_HEIGHT, Z0)

    # uref_measurement = uref_measurement_norm * ustar_ref_model

    filename = 'postProcessing/verticalProfiles/{time}/{profile}_scalar.xy'.format(
        time=latestTime,
        profile=REF_PROFILE
    )

    z, epsilon, k, p = read_epsilon_k_p(filenameS)
    # k_ref_model = k[ref_height_index]
    TIref_model = np.sqrt(2 / 3 * k) / (Uref_mag_model + 1e-30)

    # iref = np.sqrt(2 / 3 * k_ref_model) / uref_model

    U_axis = (ax1, ax3, ax5, ax7)
    for i, profile_name in enumerate(PROFILE_NAMES):
        filename = 'postProcessing/verticalProfiles/{time}/{profile}_vector.xy'.format(
            time=latestTime,
            profile=profile_name
        )
        z, U = read_U(filenameV)
        Umag = np.linalg.norm(U, axis=1)

        # measured data is already normalized by ustar and tree height
        U_measured = np.extract(
            (measurements['Mast'] == i + 1),
            measurements
        )['U_n']

        z_measured = np.extract(
            (measurements['Mast'] == i + 1),
            measurements
        )['z_n']

        U_axis[i].plot(
            U[:, 0] / ustar_ref_model,
            z / TREE_HEIGHT, label='model'
        )

        U_axis[i].plot(
            U_measured, z_measured,
            '*',
            label='measured',
        )

        U_axis[i].minorticks_on()
        U_axis[i].set_ylim((0, 5))
        U_axis[i].set_ylabel('z/h')
        U_axis[i].set_title('%s' % profile_name)
        U_axis[i].grid(visible=True, which='minor', color='black',
                       linestyle='-', linewidth=0.1)
        U_axis[i].grid(visible=True, which='major', color='black',
                       linestyle='-', linewidth=0.5)
    U_axis[i].set_xlabel('U/$u*_{ref}$')

    TI_axis = (ax2, ax4, ax6, ax8)
    for i, profile_name in enumerate(PROFILE_NAMES):
        filename = 'postProcessing/verticalProfiles/{time}/{profile}_scalar.xy'.format(
            time=latestTime,
            profile=profile_name
        )
        z, epsilon, k, p = read_epsilon_k_p(filenameS)

        filename = 'postProcessing/verticalProfiles/{time}/{profile}_vector.xy'.format(
            time=latestTime,
            profile=profile_name
        )
        z, U = read_U(filenameV)
        Umag = np.linalg.norm(U, axis=1)

        TI = np.sqrt(2 / 3 * k) / (Umag + 1e-30)

        z_n = np.extract(
            (measurements['Mast'] == i + 1),
            measurements
        )['z_n']

        sigmau = np.extract(
            (measurements['Mast'] == i + 1),
            measurements
        )['sigmau_n'] * ustar_ref_model

        sigmaw = np.extract(
            (measurements['Mast'] == i + 1),
            measurements
        )['sigmaw_n'] * ustar_ref_model

        sigmav = 2.1 * ustar_ref_model

        U_measured = np.extract(
            (measurements['Mast'] == i + 1),
            measurements)['U_n'] * ustar_ref_model

        W_measured = np.extract(
            (measurements['Mast'] == i + 1),
            measurements)['W_n'] * ustar_ref_model

        Umag_measured = np.sqrt(
            np.power(U_measured, 2) + np.power(W_measured, 2)
        )
        TI_measured = calculate_TI(sigmau, sigmav, sigmaw, Umag_measured)

        TI_axis[i].plot(
            TI,
            z / TREE_HEIGHT, label='model'
        )

        TI_axis[i].plot(
            TI_measured,
            z_n,
            '*', label='measured'
        )
        TI_axis[i].minorticks_on()
        TI_axis[i].grid(visible=True, which='minor', color='black',
                        linestyle='-', linewidth=0.1)
        TI_axis[i].grid(visible=True, which='major', color='black',
                        linestyle='-', linewidth=0.5)

        TI_axis[i].set_ylim((0, 5))
        TI_axis[i].set_xlim((0, 2.5))
        TI_axis[i].set_title('%s' % profile_name)
    TI_axis[i].set_xlabel('TI')

    ax1.annotate(
        'Wind speed', xy=(0.5, 1), xytext=(0, 20),
        xycoords='axes fraction', textcoords='offset points',
        size='large', ha='center', va='baseline'
    )
    ax2.annotate(
        'Turbulence intensity', xy=(0.5, 1), xytext=(0, 20),
        xycoords='axes fraction', textcoords='offset points',
        size='large', ha='center', va='baseline'
    )

    # plt.show()
    plt.savefig('profiles.pdf')


if __name__ == '__main__':
    main()
