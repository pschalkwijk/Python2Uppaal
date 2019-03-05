#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 11:34:40 2018

@author: ggleizer
"""

import numpy as np
import etctime as etc
import time
import pprint

'''Plant, controller and PETC implementation from Heemels et al (2013) - Periodic Event-Triggered 
Control for Linear Systems'''

# Plant
Ap = np.array([[0, 1], [-2, 3]])
Bp = np.array([0, 1])
E = np.array([1, 0])
Cp = np.eye(2)

# Controller
K = np.array([[1, -4]])
h = 0.05

# Performance requirements
rho = 0.01  # decay rate
gamma = 2   # L-2 gain

# Triggering coefficient
sigma = 0.1

# Provided Lyapunov Matrix
Pl = np.array([[1, .25], [.25, 1]])

# Initial conditions for simulation
x0 = np.array([1, -1])

plant = etc.LinearPlant(Ap, Bp, Cp, None, E)
controller = etc.LinearController(K, h)
trig = etc.TabuadaPETC(plant, controller, Pl, None, sigma, 0, None, 5)
t = time.time()
traffic = etc.TrafficModelPETC(trig)
print('Elapsed: %.2f seconds' % (time.time() - t))

costrate = {(i, k): traffic.cost[(i, k)] / (k + 1) for i, k in traffic.cost}

from abstractions import AbstractedTA
test = AbstractedTA(traffic)
pprint.pprint(test.invariants)

from NTA import SigmaNTA
sNTA = SigmaNTA(test)

from ControlLoop import ControlLoop
cl = ControlLoop(sNTA)

from Network import Network
net = Network(2,1)

from NTA import NTGA
nta = NTGA(net, cl)


with open("test.xml", 'w') as fw:
    fw.write(nta.to_xml())
print('Elapsed %.2f seconds' % (time.time() - t))
