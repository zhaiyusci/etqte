#!/usr/bin/env python3

# ==============================================================================
#       MM\ MMMMMMMM\ MMMMMMMM\ MM\  MMMMMM\  MM\ MMMMMMMM\ MMMMMMMM\ MM\      =
#      MM  |MM  _____|\__MM  __|MM |MM  __MM\ MM |\__MM  __|MM  _____|\MM\     =
#     MM  / MM |         MM |   MM |MM /  MM |MM |   MM |   MM |       \MM\    =
#    MM  /  MMMMM\       MM |   MM |MM |  MM |MM |   MM |   MMMMM\      \MM\   =
#    \MM<   MM  __|      MM |   MM |MM |  MM |MM |   MM |   MM  __|     MM  |  =
#     \MM\  MM |         MM |   MM |MM MM\MM |MM |   MM |   MM |       MM  /   =
#      \MM\ MMMMMMMM\    MM |   MM |\MMMMMM / MM |   MM |   MMMMMMMM\ MM  /    =
#       \__|\________|   \__|   \__| \___MMM\ \__|   \__|   \________|\__/     =
#                                        \___|                                 =
# ==============================================================================
# ==============================================================================
# ========== Energy    Transfer       Quantum        Time   Evolution ==========
# ==============================================================================
# =========================== Hamiltonian Extraction ===========================
# ==============================================================================
# ================================ Version 1.0.0 ===============================
# ================================== April 2023 ================================
# ==============================================================================
# ===================== J.-R. Li, Y. Zhai, Z. Qu, and H. Li ====================
# ==============================================================================
# This program performs the time evolution based on quantum mechanics.
# The detailed theory can be found in the manual.
# The advantage of this program is that it is free from complex computations,
# and utilize Chebyshev propagator to get the population evolution.
#
# This program is released under MIT licence.
#
# Copyright 2023 the authors of ETQTE
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the “Software”), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import re
import numpy as np
import os
from typing import Any
from copy import copy
from sklearn.decomposition import PCA
from tqdm import tqdm

#  np.set_printoptions(precision=5, suppress=True)

# Global settings
#  pigments = ["5620", "5621", "602", "603",
#  "604", "610", "611", "612", "613"]
pigments = ["5620", "5621",
            "602", "603", "604", "610", "611", "612", "613",
            "602", "603", "604", "610", "611", "612", "613", ]
states = [1, 1,
          1, 1, 1, 1, 1, 1, 1,
          2, 2, 2, 2, 2, 2, 2,]
diagshift = [-0.1027, -0.1486,
             -0.0954, -0.1237, -0.0066, -0.0211, -0.019, 0.0068, -0.011,
             -0.1704, -0.1291, 0.0022, -0.0564, -0.0206, -0.0144, -0.0286,]
nstate = len(pigments)


def kabsch(xe, x):
    xec = copy(xe)
    xc = copy(x)
    xec -= np.average(xec, 0)
    xc -= np.average(xc, 0)
    c = xec.transpose() @ xc
    w, sigma, vh = np.linalg.svd(c)
    v = vh.transpose()
    d = np.sign(np.linalg.det(v@w.transpose()))
    u = v@np.diag([1, 1, d])@w.transpose()
    return u


#  with open(fname, 'r') as f:
    #  fulltext = f.read()

def findhamiltonian(fulltext, state1, state2):
    "Find the 2-state Frenkel exciton Hamiltonian from a Gaussian 16 EET computation."
    H = np.zeros((2, 2))
    # The following RegEx is based on output log file of Gaussian 16
    # A sample of it is like
    '''
 Frag=  2 State=  1 (w=  2.0972 eV) <=> Frag=  1 State=  1 (w=  2.0721 eV)

   delta-w                   =  0.025170439 eV
   Coulomb                   = -0.000608254 eV
   Exact-exchange            =  0.000000000 eV
   Exchange-correlation      =  0.000000000 eV
   w-avg*Overlap             = -0.000006478 eV (w-avg=  2.0846 eV, Ovlp=-0.310732D-05)
   Total coupling            = -0.000614731 eV
    '''
    #                                       1                                        2                               3
    #                                     H(1,1)                                  H(0,0)                          H(0,1)
    pattern = r'Frag=  2 State=  {} \(w= *(\S+?) eV\) <=> Frag=  1 State=  {} \(w= *(\S+?) eV\).*?Total coupling *= *(\S+?) eV'.format(
        state2, state1)
    resH = re.search(pattern, fulltext, re.MULTILINE | re.DOTALL)
    # Because the arbitary phase f the wave functions,
    # we need the transition dipole moment to check the sign of off-diagonal elements
    # A sample of it is like
    '''
 EET: doing calculation for fragment   1 using the basis set of the fragment.
    '''
    # and the information is in
    '''
 Ground to excited state transition electric dipole moments (Au):
       state          X           Y           Z        Dip. S.      Osc.
         1         4.6256      4.5599     -2.0110     46.2331      2.7624
         2        -0.8366     -1.0398      0.3510      1.9043      0.1368
         3        -0.7886      0.3645      0.4127      0.9250      0.0802
 Ground to excited state transition velocity dipole moments (Au):
        '''

    if (resH != None):
        H[0, 0], H[1, 1], H[0, 1] = float(
            resH[2]), float(resH[1]), float(resH[3])
    else:
        H[0, 0], H[1, 1], H[0, 1] = 0.0, 0.0, 0.0

    return H


def finddipole(fulltext, state1, state2):
    "Return the transition dipole moments of the specific states in standard orientation."

    pattern1 = r'''EET: doing calculation for fragment   1 using the basis set of the fragment.*?
 Ground to excited state transition electric dipole moments \(Au\):
       state          X           Y           Z        Dip\. S\.      Osc\.
 (.*?)
 Ground to excited state transition velocity dipole moments \(Au\):'''
    res1 = re.search(
        pattern1, fulltext, re.MULTILINE | re.DOTALL)
    if res1 != None:
        dipolelines1 = res1[1]
    else:
        dipolelines1 = ''

    pattern2 = r'''EET: doing calculation for fragment   2 using the basis set of the fragment.*?
 Ground to excited state transition electric dipole moments \(Au\):
       state          X           Y           Z        Dip\. S\.      Osc\.
 (.*?)
 Ground to excited state transition velocity dipole moments \(Au\):'''
    res2 = re.search(
        pattern2, fulltext, re.MULTILINE | re.DOTALL)
    if res2 != None:
        dipolelines2 = res2[1]
    else:
        dipolelines2 = ''

    dipole1 = np.array(
        list(map(float, dipolelines1.split('\n')[state1 - 1].split()[1:4])))
    dipole2 = np.array(
        list(map(float, dipolelines2.split('\n')[state2 - 1].split()[1:4])))
    return dipole1, dipole2


def findusergeo(fulltext):
    "Return the geometry in the input orientation."
    pattern_user_geo1 = r''' Charge =  0 Multiplicity = 1 in fragment      2\.
(.*?)
[^\n]*?Fragment=2'''
    usergeo1txt = re.search(pattern_user_geo1, fulltext,
                            re.MULTILINE | re.DOTALL)
    if usergeo1txt:
        usergeo1txt = usergeo1txt[1].split('\n')
    else:
        usergeo1txt = []

    pattern_user_geo2 = r'''Fragment=1\) [^\n]*
[^\n]*?Fragment=2\) .*?
 Stoichiometry'''
    usergeo2txt = re.search(pattern_user_geo2, fulltext,
                            re.MULTILINE | re.DOTALL)
    if usergeo2txt:
        usergeo2txt = usergeo2txt[0].split('\n')[1:-2]
    else:
        usergeo2txt = []

    #  print(usergeo1txt)
    #  print(usergeo2txt)

    usergeo1 = np.zeros((len(usergeo1txt), 3))
    usergeo2 = np.zeros((len(usergeo2txt), 3))
    for i in range(len(usergeo1txt)):
        words = usergeo1txt[i].split()
        if words:
            usergeo1[i, 0], usergeo1[i, 1], usergeo1[i, 2] = float(
                words[1]), float(words[2]), float(words[3])
    for i in range(len(usergeo2txt)):
        words = usergeo2txt[i].split()
        if words:
            usergeo2[i, 0], usergeo2[i, 1], usergeo2[i, 2] = float(
                words[1]), float(words[2]), float(words[3])
    return usergeo1, usergeo2

    #  return u1.transpose()@dipole1, u2.transpose()@dipole2


def findstdgeo(fulltext):
    """Return the 'standard orientation'. 
    Note that we do not split the two fragment for tech reason.
    """

    pattern_std_geo = r'''
                         Standard orientation:                         
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates \(Angstroms\)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
 (.*)
 ---------------------------------------------------------------------
 '''
    std_geo_lines = re.search(
        pattern_std_geo, fulltext, re.MULTILINE | re.DOTALL)
    if std_geo_lines:
        std_geo_lines = std_geo_lines[1].split('\n')
    else:
        std_geo_lines = []

    stdgeo = np.zeros((len(std_geo_lines), 3))

    if std_geo_lines:
        for i in range(len(std_geo_lines)):
            words = std_geo_lines[i].split()
            stdgeo[i, 0], stdgeo[i, 1], stdgeo[i, 2] = float(
                words[3]), float(words[4]), float(words[5])

    return stdgeo


def cosanglebetween(v0, v1):
    value = v0 @ v1 / np.linalg.norm(v0) / np.linalg.norm(v1)
    return value


def hamiltoniantwobytwo(prefix, i, j):
    fname = f'{prefix}{pigments[i]}-{pigments[j]}-eet.log'
    try:
        with open(fname, "r") as f:
            fulltext = f.read()
            HH = findhamiltonian(fulltext, states[i], states[j])
            mui, muj = finddipole(fulltext, states[i], states[j])
            stdgeo = findstdgeo(fulltext)
            usergeoi, usergeoj = findusergeo(fulltext)

        ui = kabsch(stdgeo[0:usergeoi.shape[0], :], usergeoi)
        uj = kabsch(stdgeo[usergeoi.shape[0]:, :], usergeoj)

        mui = ui@mui
        muj = uj@muj

    except FileNotFoundError:
        HH, muj, mui, usergeoj, usergeoi = hamiltoniantwobytwo(prefix, j, i)
        HH[0, 0], HH[1, 1] = HH[1, 1], HH[0, 0]

    return HH, mui, muj, usergeoi, usergeoj


def extracthamiltonian(prefix, diagshift=[]):
    "Find the N-state Frenkel exciton Hamiltonian from 2-state Frenkel Hamiltonians."

    H = np.zeros((nstate, nstate))
    mu = np.zeros((nstate, 3))
    mu.fill(np.nan)
    geo: list[Any] = [None for _ in range(nstate)]

    for i in range(nstate):
        for j in range(i+1, nstate):
            if pigments[i] != pigments[j]:
                HH, mui, muj, usergeoi, usergeoj = hamiltoniantwobytwo(
                    prefix, i, j)
                if geo[i] is None:
                    geo[i] = copy(usergeoi)
                if geo[j] is None:
                    geo[j] = copy(usergeoj)

                if np.isnan(mu[i, 0]):
                    mu[i, :] = mui[:]
                else:
                    cc = cosanglebetween(mu[i, :], mui)
                    if -0.95 < cc < 0.95:  # Should be +/- 1
                        print("*** cosine of the two trans. dip. is", cc)
                    if cc < 0:
                        HH[0, 1] *= -1.0

                if np.isnan(mu[j, 0]):
                    mu[j, :] = muj[:]
                else:
                    cc = cosanglebetween(mu[j, :], muj)
                    if -0.95 < cc < 0.95:  # Should be +/- 1
                        print("*** cosine of the two trans. dip. is", cc)
                    if cc < 0:
                        HH[0, 1] *= -1.0

                #  if H[i,i] != 0.0:
                    #  print("Delta H = ", H[i,i] - HH[0,0])
                H[np.ix_([i, j], [i, j])] = copy(HH)

    if len(diagshift) != 0:
        H += np.diag(diagshift)

    return H, mu, geo


def phase1():
    """
    This function serves as the phase correct for each isolated system.
    """
    try:
        os.mkdir("heff1")
    except OSError as error:
        print(error)

    for i in range(nstate):
        with open(f"{i}.dipole.tmp", 'w') as f:
            f.write('')
        with open(f"geo{i}.xyz", 'w') as f:
            f.write('')

    for s in tqdm(range(1, 1001)):
        H, mu, geo = extracthamiltonian(f"../rawdata/{s}/", diagshift)

        np.savetxt(f"heff1/{s}.mat", np.triu(H), "%.9f")

        for i in range(nstate):
            with open(f"{i}.dipole.tmp", 'a') as f:
                f.write(f"{mu[i][0]} {mu[i][1]} {mu[i][2]}\n")
            with open(f"geo{i}.xyz", "a") as f:
                f.write(f"{geo[i].shape[0]}\n\n")
                for atom in range(geo[i].shape[0]):
                    f.write(
                        f"C   {geo[i][atom, 0]} {geo[i][atom, 1]} {geo[i][atom, 2]}\n")

    print("phase1() is done.")
    return


def phase2():
    """
    This function tries align the trans. dip. of all systems in the ensemble.
    """

    try:
        os.mkdir("heff2")
    except OSError as error:
        print(error)

    dips = [np.empty((1001, 3)) for _ in range(nstate)]
    refdip = [np.empty(3) for _ in range(nstate)]
    for i in range(nstate):
        pca = PCA(n_components=3)
        dips[i] = np.loadtxt(f"{i}.dipole.tmp")
        pca.fit(np.concatenate((dips[i], -dips[i])))
        refdip[i] = copy(pca.components_[0])
    with open("refdip", "w") as f:
        for dip in refdip:
            f.write(f"{dip[0]*5} {dip[1]*5} {dip[2]*5} \n\n ")

    with open("total.stath", "w") as f:
        f.write('')

    for s in tqdm(range(1, 1001)):
        H = np.loadtxt(f"heff1/{s}.mat")
        thedip = [dips[i][s-1] for i in range(nstate)]
        cc = np.array([cosanglebetween(refdip[i], thedip[i])
                      for i in range(nstate)])
        cc[cc > 0] = 1
        cc[cc < 0] = -1
        H = np.diag(cc) @ H @ np.diag(cc)
        np.savetxt(f"heff2/{s}.mat", np.triu(H), "%.9f")
        with open("total.stath", "a") as f:
            for i in range(nstate):
                for j in range(i, nstate):
                    f.write(f"{H[i,j]} ")
            f.write("\n")
    print("phase2() is done.")
    return


if __name__ == '__main__':
    phase1()
    phase2()
    # test
