import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from dataPlotterObserver import dataPlotterObserver
from blockBeamDynamics import blockBeamDynamics
from ctrlPID import ctrlPID
from ctrlObserver14 import ctrlObserver

# instantiate blockBeam, controller, and reference classes
blockBeam = blockBeamDynamics(alpha=0.)
controller = ctrlObserver()
reference = signalGenerator(amplitude=.125, frequency=0.05, y_offset=.25)
disturbance = signalGenerator(amplitude=.25)
noise_z = signalGenerator(amplitude=.01)
noise_th = signalGenerator(amplitude=.01)


# instantiate the simulation plots and animation
dataPlot = dataPlotter()
dataPlotterObserver = dataPlotterObserver()
animation = blockbeamAnimation()

t = P.t_start  # time starts at t_start
y = blockBeam.h()  # output of system at start of simulation
while t < P.t_end:  # main simulation loop
    # Get referenced inputs from signal generators
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot

    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        r = reference.square(t)
        d = 0.0#disturbance.step(t)  # input disturbance
        n = 0.#np.array([[noise_z.random(t)], [noise_th.random(t)]])  # simulate sensor noise
        x = blockBeam.state
        u, xhat, dhat = controller.update(r, y + n)  # update controller
        y = blockBeam.update(u + d)  # propagate system
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(blockBeam.state)
    dataPlot.update(t, r, blockBeam.state, u)
    dataPlotterObserver.update(t, blockBeam.state, xhat, d, dhat)

    # the pause causes the figure to be displayed for simulation
    #plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
#dataPlot.write_data_file()
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
