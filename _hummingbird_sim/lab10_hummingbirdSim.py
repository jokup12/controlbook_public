import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import HummingbirdDynamics
from controller_10 import ctrlPID

# instantiate pendulum, controller, and reference classes
hummingbird = HummingbirdDynamics(alpha=0.)
controller = ctrlPID()
phi_ref = SignalGenerator(amplitude=15.*np.pi/180., frequency=0.02)
theta_ref = SignalGenerator(amplitude=15.*np.pi/180., frequency=0.05)
psi_ref = SignalGenerator(amplitude=30.*np.pi/180., frequency=0.02)
disturbance = SignalGenerator(amplitude=.02)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()

t = P.t_start  # time starts at t_start
y = hummingbird.h()
while t < P.t_end:  # main simulation loop

    # Propagate dynamics at rate Ts
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = np.array([[0], [theta_ref.square(t)], [psi_ref.square(t)]])
        u, y_ref, v = controller.update(r, y)
        force = v[0][0]
        torque = v[1][0]

        d = force*.02#disturbance.step(t)
        y = hummingbird.update(u + d)  # Propagate the dynamics
        t = t + P.Ts  # advance time by Ts



        y_ref[2][0] = psi_ref.square(t)

    # update animation and data plots at rate t_plot
    animation.update(t, hummingbird.state)
    dataPlot.update(t, hummingbird.state, y_ref, force, torque)

    # the pause causes figure to be displayed during simulation
    #plt.pause(0.00001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
