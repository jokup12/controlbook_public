import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter

# instantiate reference input classes
phi_ref = SignalGenerator(amplitude=1, frequency=0.3)
theta_ref = SignalGenerator(amplitude=1.5, frequency=0.1)
psi_ref = SignalGenerator(amplitude=2, frequency=0.15)
force_ref = SignalGenerator(amplitude=10, frequency=.4, y_offset=20)
torque_ref = SignalGenerator(amplitude=.5, frequency=.23, y_offset=1)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    phi = phi_ref.sin(t)            #roll, positive spin is forward
    theta = theta_ref.sin(t)      #pitch, positive spin is starboard
    psi = psi_ref.sin(t)          #yaw, positive spin is down
    # update animation
    state = np.array([[phi], [theta], [psi], [0.0], [0.0], [0.0]])
    ref = np.array([[0], [0], [0]])
    force = force_ref.square(t)
    torque = torque_ref.square(t)
    animation.update(t, state)
    dataPlot.update(t, state, ref, force, torque)

    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.05)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()