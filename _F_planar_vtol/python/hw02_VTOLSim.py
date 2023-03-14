import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter

# instantiate reference input classes
position = signalGenerator(amplitude=4, frequency=0.1, y_offset=5)
altitude = signalGenerator(amplitude=2.0, frequency=0.1, y_offset=2)
angle = signalGenerator(amplitude=np.pi/8, frequency=.5)
force = signalGenerator(amplitude=np.pi/8, frequency=.5)
torque = signalGenerator(amplitude=np.pi/8, frequency=.5)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    z = position.sin(t)
    h = altitude.square(t)
    theta = angle.sin(t)
    F = force.square(t)
    T = torque.square(t)
    # update animation
    state = np.array([[z], [h], [theta]])  # state is made of theta, and theta_dot
    animation.update(state)
    dataPlot.update(t, state, z, h, F, T)
    # advance time by t_plot
    t = t + P.t_plot
    plt.pause(0.01)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
