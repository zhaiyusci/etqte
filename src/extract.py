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


def findhamiltonian(fname, state1, state2):
    "Find the 2-state Frenkel exciton Hamiltonian from a Gaussian 16 EET computation."
    with open(fname, 'r') as f:
        fulltext = f.read()
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
        res = re.search(pattern, fulltext, re.MULTILINE | re.DOTALL)
    if (res != None):
        return float(res[2]), float(res[1]), float(res[3])
    else:
        return 0.0, 0.0, 0.0


def findtransitiondipole(fname, state1, state2):
    "Find the transition dipole moment from ground state to selected states from a Gaussian 16 EET computation."
    with open(fname, 'r') as f:
        fulltext = f.read()
        # The following RegEx is based on output log file of Gaussian 16
        # A sample of it is led by
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

        pattern1 = r'''EET: doing calculation for fragment   1 using the basis set of the fragment.*?
 Ground to excited state transition electric dipole moments \(Au\):
       state          X           Y           Z        Dip\. S\.      Osc\.
(.*?)
 Ground to excited state transition velocity dipole moments \(Au\):'''
        res = re.search(
            pattern1, fulltext, re.MULTILINE | re.DOTALL)
        if res != None:
            dipolelines1 = res[1]
        else:
            dipolelines1 = ''

        pattern2 = r'''EET: doing calculation for fragment   2 using the basis set of the fragment.*?
 Ground to excited state transition electric dipole moments \(Au\):
       state          X           Y           Z        Dip\. S\.      Osc\.
(.*?)
 Ground to excited state transition velocity dipole moments \(Au\):'''
        res = re.search(
            pattern2, fulltext, re.MULTILINE | re.DOTALL)
        if res != None:
            dipolelines2 = res[1]
        else:
            dipolelines2 = ''

        dipole1 = np.array(
            list(map(float, dipolelines1.split('\n')[state1 - 1].split()[1:4])))
        dipole2 = np.array(
            list(map(float, dipolelines2.split('\n')[state2 - 1].split()[1:4])))
    return dipole1, dipole2


def extracthamiltonian(prefix, pigments, states, diagshift=[]):
    "Find the N-state Frenkel exciton Hamiltonian from 2-state Frenkel Hamiltonians."
    if len(pigments) == len(states):
        nstate = len(pigments)
        H = np.empty((nstate, nstate))
    else:
        exit(-1)

    for i in range(9):
        for j in range(i+1, 9):
            fname = f'{prefix}{pigments[i]}-{pigments[j]}-eet.log'
            #  print(fname)

            # We only need the upper triangle matrix.
            H[i, i], H[j, j], H[i, j] = findhamiltonian(
                fname, states[i], states[j])
    if len(diagshift) != 0:
        H += np.diag(diagshift)

    return H


def main():
    pigments = ["5620", "5621", "602", "603",
                "604", "610", "611", "612", "613"]
    states = [1, 1, 2, 2, 2, 2, 2, 2, 2]
    diagshift = [-0.1027, -0.1486, -0.0954, -0.1237, -0.0066, -0.0211, -0.019, 0.0068, -0.011]
    try:
        os.mkdir("heff")
    except OSError as error:
        print(error)

    for i in range(1, 4):
        #  print(i)
        H = extracthamiltonian(f"{i}/", pigments, states, diagshift)
        with open("H00.stat", 'a') as f:
            f.write(f"{H[0,0]}\n")
        with open("H11.stat", 'a') as f:
            f.write(f"{H[1,1]}\n")
        with open("H01.stat", 'a') as f:
            f.write(f"{H[0,1]}\n")
        with open(f"heff/{i}.mat", 'w') as f:
            for i in range(len(pigments)):
                for j in range(i, len(pigments)):
                    f.write(
                        f'{pigments[i]:10s}{pigments[j]:10s}{H[i,j]:20.10f}\n')

    return


if __name__ == '__main__':
    main()
    # test
    states = [1, 1, 2, 2, 2, 2, 2, 2, 2]
    for i in range(1, 4):
        print(findtransitiondipole(
            f'{i}/5620-5621-eet.log', states[0], states[1]))
