import matplotlib.pyplot as plt
import VTOLParam as P
import numpy as np
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics

# instantiate VTOL, controller, and reference classes
VTOL = VTOLDynamics()
# second_VTOL = VTOLDynamics()
reference = signalGenerator(amplitude=10, frequency=1)
force = signalGenerator(amplitude=10, frequency=1)
# torque2 = signalGenerator(amplitude=0.1, frequency=0.01)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
# dataPlot2 = dataPlotter()
animation = VTOLAnimation()
# animation2 = VTOLAnimation()


t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        # Get referenced inputs from signal generators
        r = reference.sin(t)
        #u[0][0] = force.sin(t)
        #u[1][0] = force.sin(t)
        u = np.array([[force.sin(t)], [force.sin(t)]])
        f = force.square(t)
        tor = force.sawtooth(t)
        # u2 = torque2.sin(t)
        y = VTOL.update(u)  # Propagate the dynamics
        # y2 = second_VTOL.update(u2)
        t = t + P.Ts  # advance time by Ts
    # update animation and data plots
    animation.update(VTOL.state)
    # animation2.update(second_VTOL.state)
    dataPlot.update(t, VTOL.state, r, r, f, tor)
    # dataPlot2.update(t, r, second_VTOL.state, u2)

    # the pause causes the figure to be displayed during the
    # simulation
    plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
