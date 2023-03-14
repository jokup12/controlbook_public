import matplotlib.pyplot as plt
import VTOLParam as P
import numpy as np
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics
from ctrlPID import ctrlPID

# instantiate VTOL, controller, and reference classes
VTOL = VTOLDynamics(alpha=.2)
controller = ctrlPID()
reference = signalGenerator(amplitude=2.5, frequency=.08, y_offset=3)
height = signalGenerator(amplitude=1, frequency=.03, y_offset=2)
force = signalGenerator(amplitude=10, frequency=1)
# torque2 = signalGenerator(amplitude=0.1, frequency=0.01)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
# dataPlot2 = dataPlotter()
animation = VTOLAnimation()
# animation2 = VTOLAnimation()


t = P.t_start  # time starts at t_start
y = VTOL.h()
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        # Get referenced inputs from signal generators
        r = reference.square(t)
        h = height.square(t)
        x = VTOL.state
        u, v = controller.update(r, h, x)  # update controller
        #print(u)
        f = u[0][0]
        tor = u[1][0]



        #f = np.array([[fr], [fl]])

        y = VTOL.update(v)  # Propagate the dynamics

        t = t + P.Ts  # advance time by Ts
    # update animation and data plots

    animation.update(VTOL.state)
    dataPlot.update(t, VTOL.state, r, h, f, tor)

    # the pause causes the figure to be displayed during the
    # simulation
    plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
