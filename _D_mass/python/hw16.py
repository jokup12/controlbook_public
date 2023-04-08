# Single link arm Parameter File
import massParam as P
from ctrlPID import ctrlPID
import hw15 as P15
from control import tf, bode
import matplotlib.pyplot as plt

P10 = ctrlPID()
sigma = .05

# flag to define if using dB or absolute scale for M(omega)
dB_flag = False # False

# Assign plan from previous homework solution
Plant = P15.Plant

# Compute transfer function of controller
C_pid = tf([(P10.kd+P10.kp*sigma), (P10.kp+P10.ki*sigma), P10.ki],
           [sigma, 1, 0])

if __name__ == '__main__':
    # display bode plots of transfer functions
    fig = plt.figure()
    bode([Plant, Plant*C_pid], dB=dB_flag)
    fig.suptitle('Mass Spring Damper')
    plt.legend(['P(s)', 'C(s)P(s)'])
    print('Close plot window to end program')
    plt.show()
