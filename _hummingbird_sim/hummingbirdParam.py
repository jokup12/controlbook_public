# Hummingbird Parameter File
import numpy as np
# Initial Conditions
phi0 = 0.0 * np.pi / 180  # roll angle in rads
theta0 = 0 * np.pi / 180  # pitch angle in rads
psi0 = 0.0 * np.pi / 180  # yaw angle in rads
phidot0 = 0.0              # roll rate in rads/sec
thetadot0 = 0.0         # pitch rate in rads/sec
psidot0 = 0.0              # yaw rate in rads/sec
# Physical parameters of the hummingbird known to the controller
g = 9.81
ell1 = .247
ell2 = -.039
ell3x = -.007
ell3y = -.007
ell3z = .018
ellT = .355
d = .12
m1 = .108862
J1x = .000189
J1y = .001953
J1z = .001894
m2 = .4717
J2x = .00231
J2y = .003274
J2z = .003416
m3 = .1905
J3x = .0002222
J3y = .0001956
J3z = .000027
km = g * (m1 * ell1 + m2 * ell2) / ellT  # need to find this experimentally for hardware
btheta = .001
bphi = .001
bpsi = .001


# Simulation Parameters
t_start = 0.0  # Start time of simulation
t_end = 100.0  # End time of simulation
Ts = 0.01  # sample time for simulation
t_plot = 0.1  # the plotting and animation is updated at this rate

# dirty derivative parameters
sigma = .005  # cutoff freq for dirty derivative
# beta =  # dirty derivative gain

# saturation limits
Fmax = 50.0               # Max force N
Fmin = 0.0
PWM_max = .7
PWM_min = 0.0
torque_max = 5.0                # Max torque, Nm