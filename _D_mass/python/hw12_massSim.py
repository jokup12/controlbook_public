import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics
from ctrlPID import ctrlPID
from ctrlStateFeedbackIntegrator import ctrlStateFeedback

# instantiate mass, controller, and reference classes
mass = massDynamics(alpha=0.2)
controller = ctrlStateFeedback()
reference = signalGenerator(amplitude=1, frequency=0.05)
disturbance = signalGenerator(amplitude=0.25)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()

t = P.t_start  # time starts at t_start
y = mass.h()  # output of system at start of simulation
while t < P.t_end:  # main simulation loop
    # Get referenced inputs from signal generators
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot

    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        r = reference.square(t)
        d = disturbance.step(t)  # input disturbance
        n = 0.0  #noise.random(t)  # simulate sensor noise
        x = mass.state
        u = controller.update(r, mass.state)  # update controller
        y = mass.update(u + d)  # propagate system
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(mass.state)
    dataPlot.update(t, r, mass.state, u)

    # the pause causes the figure to be displayed for simulation
    plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
dataPlot.write_data_file()
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
