#!/usr/bin/python3
# --------------------------------------------------------------------------- #
# |       o        |
# |    o     o     |  HELYX (R) : Open-source CFD for Enterprise
# |   o   O   o    |  Version : 4.4.0
# |    o     o     |  ENGYS Ltd. <http://engys.com/>
# |       o        |
# --------------------------------------------------------------------------- #
# Description
#     Plot validation graphs for different windkessel initial inputs.
#
# --------------------------------------------------------------------------- #

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def windkessel(Q, t, Rp, C, Rd, p0=0):
    for i, ti in enumerate(t):
        if i == 0:
            p[i] = p0
        else:
            dt = t[i] - t[i-1]
            p[i] = (p[i-1] + (dt/C + Rp*dt/(C*Rd) + Rp)*Q[i] - Rp*Q[i-1])/(1.0 + dt/(C*Rd))

    return p

def plotVerification(index, t, p, Q, fileName):
    CSV = np.genfromtxt(fileName, delimiter=" ", skip_header=1)
    tf = CSV[:, 0]
    Qf = CSV[:, 1]
    pf = CSV[:, 8]

    rho = 998.2

    plt.figure()
    plt.clf()

    colors = ["orange", "green", "red"]

    ax0 = plt.subplot(2, 1, 1)
    ax0.title.set_text("Input " + str(index+1))
    ax0.plot(t, p, label="Analytical")
    ax0.plot(tf[::10],
             pf[::10],
             label="HELYX",
             marker="o",
             linewidth=0,
             markersize=4,
             color=colors[index])

    ax0.set_xlim([0., 0.04])
    ax0.grid()
    ax0.legend()
    ax0.set_ylabel("p (N)")

    ax1 = plt.subplot(2, 1, 2, sharex=ax0)
    ax1.plot(t, Q)
    ax1.plot(tf[::10],
             Qf[::10]/rho,
             label="HELYX",
             marker="o",
             linewidth=0,
             markersize=4,
             color=colors[index])

    ax1.grid()
    ax1.set_xlabel("t (s)")
    ax1.set_ylabel("Q (m$^3$/s)")

    plt.setp(ax0.get_xticklabels(), visible=False)

    pp.savefig()

# Variables array
t = np.arange(0.0, 0.04, 0.0001)

Qi = interp1d(np.array([0.0, 0.01, 0.02, 0.03, 0.04]),
              np.array([1.0,  2.0, -1.0, 1.0, 1.0]))

Q = Qi(t)

p = np.zeros_like(t)

# Plot file
pp = PdfPages('windkesselValidation.pdf')

# --------------------------------------------------------------------------- #

fileName = "../postProcessing/SR_outlet.input1/0/surfaceReport.dat"

p0 = 0.0
Rp = 0
C  = 1e-5
Rd = 1e3

p = windkessel(Q, t, Rp, C, Rd, p0)

plotVerification(0, t, p, Q, fileName)

del fileName, p0, Rp, C, Rd

# --------------------------------------------------------------------------- #

fileName = "../postProcessing/SR_outlet.input2/0/surfaceReport.dat"

p0 = 500.0
Rp = 5e2
C  = 2e-5
Rd = 5e2

p = windkessel(Q, t, Rp, C, Rd, p0)

plotVerification(1, t, p, Q, fileName)

del fileName, p0, Rp, C, Rd

# --------------------------------------------------------------------------- #

Qi = interp1d(np.array([0.0, 0.00999, 0.01001, 0.0199, 0.02001, 0.02999, 0.03001, 0.04]),
              np.array([1.0,     1.0,     2.0,    2.0,    -1.0,    -1.0,     1.0, 1.0]))
Q = Qi(t)

fileName = "../postProcessing/SR_outlet.input3/0/surfaceReport.dat"

p0 = 0.0
Rp = 0
C  = 5e-6
Rd = 1e3

p = windkessel(Q, t, Rp, C, Rd, p0)

plotVerification(2, t, p, Q, fileName)

del fileName, p0, Rp, C, Rd

pp.close()

# --------------------------------------------------------------------------- #
