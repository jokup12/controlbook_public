import matplotlib.pyplot as plt
import VTOLParam as P
import numpy as np
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics
from ctrlPID import ctrlPID
from ctrlStateFeedbackIntegrator import ctrlStateFeedback

# instantiate VTOL, controller, and reference classes
VTOL = VTOLDynamics(alpha=0.2)
controller = ctrlStateFeedback()
reference = signalGenerator(amplitude=4.0, frequency=.02, y_offset=5.0)
height = signalGenerator(amplitude=3.0, frequency=.03, y_offset=5.0)
#force = signalGenerator(amplitude=10, frequency=1)
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
        ref = np.array([[r], [h]])
        d = np.array([[0.1], [0.0]])
        n = np.array([[0.0], [0.0], [0.0]])
        x = VTOL.state
        if(VTOL.state[1]<0):
            VTOL.state[1] = 0
        u = controller.update(ref, x)  # update controller
        y = VTOL.update(P.mixing @ (u+d))  # Propagate the dynamics

        #print(u)
        f = u[0][0]
        tor = u[1][0]


        t = t + P.Ts  # advance time by Ts
    # update animation and data plots

    animation.update(VTOL.state, r)
    dataPlot.update(t, VTOL.state, r, h, f, tor)

    # the pause causes the figure to be displayed during the
    # simulation
    #plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
