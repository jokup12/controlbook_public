# Single link arm Parameter File
import massParam as P
from ctrlPID import ctrlPID
import hw15 as P15
from control import tf, bode
import matplotlib.pyplot as plt
import numpy as np

P10 = ctrlPID()
sigma = .05

# flag to define if using dB or absolute scale for M(omega)
dB_flag = P15.dB_flag # False

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

    omegas = [.1, 100.0]
    mag_CP, phase, omegas = bode(Plant*C_pid, plot = False, omega = omegas)
    mag_P, phase, omegas = bode(Plant, plot=False, omega = omegas)
    plt.figure()
    xfer_func = (Plant*C_pid)._repr_latex_()
    plt.text(0.1, 0.5, '$%s$'%xfer_func[2:-2], fontsize='xx-large')
    plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False,
                    right=False, labelbottom=False, labelleft=False)

    print("\n\nvalues for metric calculation for parts b and c:")

    if dB_flag == False:
        print("gamma_d = ", mag_P[0]/mag_CP[0])
        print("gamma_n = ", mag_CP[1])

    elif dB_flag == True:
        mag_P_db = 20.0*np.log10(mag_P)
        mag_CP_dB = 20.0*np.log10(mag_CP)
        print("gamma_d = ", (mag_P[0]/mag_CP[0])/20.0)
        print("gamma_n = ", mag_CP[1]/20.0)





    print('Close plot window to end program')
    plt.show()
