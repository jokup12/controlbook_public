import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter

# instantiate reference input classes
reference = signalGenerator(amplitude=0.5, frequency=0.1)
blockPos = signalGenerator(amplitude=.05, frequency=0.5, y_offset=.2)
beamAngle = signalGenerator(amplitude=np.pi/8, frequency=.1)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    r = reference.square(t)
    z = blockPos.sin(t)
    theta = beamAngle.square(t)
    # update animation
    state = np.array([[z], [theta]])  # state is made of theta, and theta_dot
    animation.update(state)
    dataPlot.update(t, r, state, theta)
    # advance time by t_plot
    t = t + P.t_plot
    plt.pause(0.01)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
