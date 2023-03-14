import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter

# instantiate reference input classes
reference = signalGenerator(amplitude=0.5, frequency=0.1)
mass = signalGenerator(amplitude=1, frequency=0.5, y_offset=.2)
tauRef = signalGenerator(amplitude=5, frequency=.5)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    r = reference.square(t)
    massPos = mass.sin(t)
    tau = tauRef.sawtooth(t)
    # update animation
    state = np.array([[massPos], [0.0]])  # state is made of theta, and theta_dot
    animation.update(state)
    dataPlot.update(t, r, state, tau)
    # advance time by t_plot
    t = t + P.t_plot
    plt.pause(0.01)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
