#! /usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
import importlib
#plt.ion()
plt.ioff()

import grlp

#%load_ext autoreload
#%autoreload

importlib.reload(grlp)
#del net

dt = 3.15E7*10
_B = 100 # uniform

nseg = 5
numel = [5, 7, 4, 8, 6]

upstream_segment_IDs = [[], [], [0,1], [], [2,3]]
downstream_segment_IDs = [[2], [2], [4], [4], []]

z = []
Q_in_list = [5., 10., 15., 10., 25.]
Q = []
B = []
print( "" )
print( "**************" )
print( "" )
for i in range(nseg):
    # Start z as all zeros
    z.append( np.zeros( numel[i] ) )
    # Start Q as constant within each segment
    Q.append( Q_in_list[i] * np.ones( numel[i] ) )
    # Uniform valley width: Simpler for example
    B.append( _B * np.ones( numel[i] ) )
    print( numel[i] )
print( "" )
print( "**************" )
print( "" )

# Custom for just this test network
x = [
      1000 * np.array([2, 4, 6.5, 9, 10]),
      1000 * np.array([0, 1, 2, 3, 6, 8, 10.5]),
      1000 * np.array([12, 15, 18, 20]),
      1000 * np.array([2, 6, 8, 12, 14, 16, 18, 20]),
      1000 * np.array([23, 24, 27, 29, 29.5, 30])
    ]

# Base level
x_bl = 1000*32
z_bl = 0

# Upstream boundary condition: 1.5% grade
S0 = [0.03, 0.015, 0.01]

# Instantiate network object
net = grlp.Network()
self = net # For debugging, etc.

# Initialize network by passing x, z, etc. information to model
net.initialize(
                config_file = None,
                x_bl = x_bl,
                z_bl = z_bl,
                S0 = S0,
                Q_s_0 = None,
                upstream_segment_IDs = upstream_segment_IDs,
                downstream_segment_IDs = downstream_segment_IDs,
                x = x,
                z = z,
                Q = Q,
                B = B,
                overwrite=False
                )

# Should do this above
net.set_niter(1)
net.get_z_lengths()

# For testing
net.evolve_threshold_width_river_network(nt=1, dt=dt)

#"""
# For plotting
# WHEN RUN FOR NT=10, GET BACKWARDS SLOPE ON TRIBUTARY
# THIS IS WHERE WE NEED TO ADD IN CLOSED BASINS AS ANOTHER SEGMENT TYPE
net.evolve_threshold_width_river_network(nt=500, dt=dt)

for lp in net.list_of_LongProfile_objects:
    for _id in lp.downstream_segment_IDs:
        dsseg = net.list_of_LongProfile_objects[_id]
        _xjoin = [lp.x[-1], dsseg.x[0]]
        _zjoin = [lp.z[-1], dsseg.z[0]]
        plt.plot(_xjoin, _zjoin, 'k-', linewidth=4, alpha=.5)
    if len(lp.downstream_segment_IDs) == 0:
        plt.plot(lp.x_ext[1:], lp.z_ext[1:], '-', linewidth=4, alpha=.5)#, label=lp.)
    else:
        plt.plot(lp.x, lp.z, '-', linewidth=4, alpha=.5)#, label=lp.)
plt.tight_layout()
plt.show()
#"""
