import matplotlib.pyplot as plt
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics

# instantiate mass, controller, and reference classes
mass = massDynamics()
# second_mass = massDynamics()
reference = signalGenerator(amplitude=1, frequency=.2)
force = signalGenerator(amplitude=10, frequency=1)
# torque2 = signalGenerator(amplitude=0.1, frequency=0.01)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
# dataPlot2 = dataPlotter()
animation = massAnimation()
# animation2 = massAnimation()

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
        y = mass.update(u)  # Propagate the dynamics
        # y2 = second_mass.update(u2)
        t = t + P.Ts  # advance time by Ts
    # update animation and data plots
    animation.update(mass.state)
    # animation2.update(second_mass.state)
    dataPlot.update(t, r, mass.state, u)
    # dataPlot2.update(t, r, second_mass.state, u2)

    # the pause causes the figure to be displayed during the
    # simulation
    plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
