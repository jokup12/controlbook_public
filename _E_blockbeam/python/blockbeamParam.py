# Ball on Beam Parameter File
import numpy as np
# import control as cnt

# Physical parameters of the  blockbeam known to the controller
m1 = .35  # Mass of the block, kg
m2 = 2.0  # mass of beam, kg
length = .5  # length of beam, m
g = 9.8  # gravity at sea level, m/s^2

# parameters for animation
width = 0.05  # width of block
height = width*0.25  # height of block

# Initial Conditions
z0 = .25 # initial block position,m
theta0 = 0*np.pi/180  # initial beam angle,rads
zdot0 = 0.0  # initial speed of block along beam, m/s
thetadot0 = 0.0 # initial angular speed of the beam,rads/s

# Simulation Parameters
t_start = 0  # Start time of simulation
t_end = 300  # End time of simulation
Ts = .01  # sample time for simulation
t_plot = .1  # the plotting and animation is updated at this rate

# saturation limits
Fmax = 150000.0  # Max Force, N

# dirty derivative parameters
sigma = .05  # cutoff freq for dirty derivative
# beta =  # dirty derivative gain

# equilibrium force when block is in center of beam
ze = length/2
Fe = (m1*g*ze)/length + m2*g/2
