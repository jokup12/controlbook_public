import matplotlib.pyplot as plt
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockBeamDynamics import blockBeamDynamics

# instantiate blockBeam, controller, and reference classes
blockBeam = blockBeamDynamics()
# second_blockBeam = blockBeamDynamics()
reference = signalGenerator(amplitude=.5, frequency=.02)
force = signalGenerator(amplitude=.5, frequency=1, y_offset=11.5)
# torque2 = signalGenerator(amplitude=0.1, frequency=0.01)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
# dataPlot2 = dataPlotter()
animation = blockBeamAnimation()
# animation2 = blockBeamAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        # Get referenced inputs from signal generators
        r = reference.square(t)
        u = force.sin(t)
        # u2 = torque2.sin(t)
        y = blockBeam.update(u)  # Propagate the dynamics
        # y2 = second_blockBeam.update(u2)
        t = t + P.Ts  # advance time by Ts
    # update animation and data plots
    animation.update(blockBeam.state)
    # animation2.update(second_blockBeam.state)
    dataPlot.update(t, r, blockBeam.state, u)
    # dataPlot2.update(t, r, second_blockBeam.state, u2)

    # the pause causes the figure to be displayed during the
    # simulation
    plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
