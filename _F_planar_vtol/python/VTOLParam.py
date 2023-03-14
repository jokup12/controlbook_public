# VTOL Parameter File
import numpy as np

# Physical parameters of the  VTOL known to the controller
mc = 1.0 # kg
mr = .25  # kg
Jc = .0042  # kg m^2
d = .3  # m
mu = .1  # kg/s
g = 9.81  # m/s^2
F_wind = 0.0 # wind disturbance force is zero in initial homeworks

# parameters for animation
length = 10.0

# Initial Conditions
z0 = 0.0  # initial lateral position
h0 = 0.0  # initial altitude
theta0 = 0.0 # initial roll angle
zdot0 = 0.0  # initial lateral velocity
hdot0 = 0.0  # initial climb rate
thetadot0 = 0.0  # initial roll rate
target0 =0.0

# Simulation Parameters
t_start = 0 # Start time of simulation
t_end = 120  # End time of simulation
Ts = .01  # sample time for simulation
t_plot = .1 # the plotting and animation is updated at this rate

# saturation limits
Fmax = 1000000.0  # Max Force, N

# mixing matrix
mixing = np.linalg.inv(np.array([[1.0, 1.0], [d, -d]]))

# dirty derivative parameters
sigma = .05  # cutoff freq for dirty derivative
# beta =  # dirty derivative gain

# equilibrium force
Fe = (mc+2*mr)*g
Te = 0.0


