import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from hummingbirdDynamics import HummingbirdDynamics
from dataPlotter import DataPlotter
from ctrlPD_7 import ctrlPD

# instantiate reference input classes
HB = HummingbirdDynamics(alpha=.0)
controller = ctrlPD()
phi_ref = SignalGenerator(amplitude=.05, frequency=0.01)
theta_ref = SignalGenerator(amplitude=.3, frequency=0.01, y_offset=0)
psi_ref = SignalGenerator(amplitude=.05, frequency=0.01)
force_ref = SignalGenerator(amplitude=5., frequency=.4, y_offset=20)
torque_ref = SignalGenerator(amplitude=.5, frequency=.23, y_offset=1)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()

t = P.t_start  # time starts at t_start
y = HB.h()
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:    # set variables
        phi = phi_ref.square(t)            #roll, positive spin is forward
        theta = theta_ref.square(t)      #pitch, positive spin is starboard
        psi = psi_ref.square(t)          #yaw, positive spin is down
        # update animation
        state = HB.state    #phi, theta, psi, phidot, thetadot, psidot
        ref = np.array([[0], [theta], [0]])
        #force = force_ref.square(t)
        #torque = torque_ref.square(t)

        u, v = controller.update(ref, y)  # update controller
        force = u[0][0]
        torque = u[1][0]
        torque = 0.0

        y = HB.update(v)   #propogate dynamics

        t = t + P.Ts  # advance time by Ts


    animation.update(t, state)
    dataPlot.update(t, state, ref, force, torque)

    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.05)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()