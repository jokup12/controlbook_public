# mass-spring-damper Parameter File
import numpy as np

# Physical parameters of the arm known to the controller
m = 5.0  # mass kg
k = 3.0  # spring constant Kg/s^2
b = .5  # damping coefficient Kg/s
g = 9.81

# parameters for animation
length = 5.0
width = 1.0

# Initial Conditions
z0 = 0.0  # initial position of mass, m
zdot0 = 0.0  # initial velocity of mass m/s

# Simulation Parameters
t_start = 0 # Start time of simulation
t_end = 120  # End time of simulation
Ts = .01  # sample time for simulation
t_plot = .1 # the plotting and animation is updated at this rate

# dirty derivative parameters
sigma = .05 # cutoff freq for dirty derivative
# beta =   # dirty derivative gain

# saturation limits
F_max = 1000000  # Max force, N

