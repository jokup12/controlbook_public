import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockBeamDynamics import blockBeamDynamics
from ctrlPID import ctrlPID
from ctrlStateFeedback import ctrlStateFeedback

# instantiate blockBeam, controller, and reference classes
blockBeam = blockBeamDynamics(alpha=0.0)
controller = ctrlStateFeedback()
reference = signalGenerator(amplitude=.125, frequency=0.05, y_offset=.25)
disturbance = signalGenerator(amplitude=0.0)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockBeamAnimation()

t = P.t_start  # time starts at t_start
y = blockBeam.h()  # output of system at start of simulation
while t < P.t_end:  # main simulation loop
    # Get referenced inputs from signal generators
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot

    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        r = reference.square(t)
        d = 0#disturbance.step(t)  # input disturbance
        n = 0.0  #noise.random(t)  # simulate sensor noise
        x = blockBeam.state
        u = controller.update(r, x)  # update controller
        y = blockBeam.update(u + d)  # propagate system
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(blockBeam.state)
    dataPlot.update(t, r, blockBeam.state, u)

    # the pause causes the figure to be displayed for simulation
    plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
dataPlot.write_data_file()
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
