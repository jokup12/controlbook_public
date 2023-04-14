# Single link arm Parameter File
import VTOLParam as P
from control import tf, bode
import matplotlib.pyplot as plt
import hw15 as P15
import numpy as np
from ctrlPID import ctrlPID

P10 = ctrlPID()
sigma = .05

# flag to define if using dB or absolute scale for M(omega)
dB_flag = P15.dB_flag

# Compute plant transfer functions
P_lon = P15.P_lon
P_lat_in = P15.P_lat_in
P_lat_out = P15.P_lat_out

#controller transfer functions
C_lon = tf([(P10.kdH+P10.kpH*sigma), (P10.kpH+P10.kiH*sigma), P10.kiH],
           [sigma, 1, 0])
C_lat_in = tf([(P10.kdTheta+P10.kpTheta*sigma), P10.kpTheta],
           [sigma, 1])
C_lat_out = tf([(P10.kdZ+P10.kpZ*sigma), P10.kpZ],
           [sigma, 1])


if __name__ == '__main__':
    fig1 = plt.figure()
    bode([P_lon, P_lon*C_lon], dB=dB_flag, margins=False)
    fig1.axes[0].set_title('Longitude')

    fig2 = plt.figure()
    bode([P_lat_in, P_lat_in*C_lat_in], dB=dB_flag, margins=False)
    fig2.axes[0].set_title('Lateral In')

    fig3 = plt.figure()
    bode([P_lat_out, P_lat_out*C_lat_out], dB=dB_flag, margins=False)
    fig3.axes[0].set_title('Lateral Out')

    # if you want specific values at specific frequencies, you can
    # do the following (but the magnitudes are absolute, not dB)

    plt.figure()
    xfer_func = (P_lon*C_lon)._repr_latex_()
    plt.text(0.1, 0.5, '$%s$'%xfer_func[2:-2], fontsize='xx-large')
    plt.tick_params(axis='both', which='both', bottom=False, top=False,
                    left=False, right=False, labelbottom=False, labelleft=False)
    plt.title("xfer fucntion for C_lon(s)*P_lon(s)")


    #longitude controller
    omegas = [30.0]
    mag_CP_lon, phase, omegas = bode(P_lon*C_lon, plot=False, omega=omegas)

    if dB_flag == False:
        print("gammaN: ", mag_CP_lon[0])
    elif dB_flag == True:
        mag_CP_lonDB = 20.0*np.log10(mag_CP_lon)
        print("gammaN: ", 10**mag_CP_lon[0]/20)


    #inner lateral controller
    omegas = [2.0]
    mag_CP_lat_in, phase, omegas = bode(P_lat_in*C_lat_in, plot=False, omega=omegas)
    mag_P_lat_in, phase, omegas = bode(P_lat_in, plot=False, omega=omegas)

    if dB_flag == False:
        print("gammaD_in: ", mag_P_lat_in[0]/mag_CP_lat_in[0])
    elif dB_flag == True:
        mag_P_lat_in_DB = 20.0*np.log10(mag_P_lat_in)
        mag_CP_lat_in_DB = 20.0*np.log10(mag_CP_lat_in)
        print("gammaD_in: ", 10**(mag_P_lat_in_DB[0] - mag_CP_lat_in_DB[0]) /20)

    #outer lateral controller
    omegas = [.01, .1]
    mag_CP_lat_out, phase, omegas = bode(P_lat_in*C_lat_in, plot=False, omega=omegas)
    #mag_P_lat_in, phase, omegas = bode(P_lat_in, plot=False, omega=omegas)

    if dB_flag == False:
        print("gammaR: ", 1.0/mag_CP_lat_out[1])
        print("gammaD_out: ", 1.0/mag_CP_lat_out[0])
    elif dB_flag == True:
        mag_CP_lat_out_dB = 20.0*np.log10(mag_CP_lat_out)
        print("gammaR: ", 1.0 / 10**(mag_CP_lat_out_dB[0]) / 20)
        print("gammmaD_out: ", 1.0/ 10**(mag_CP_lat_out_dB[0] / 20))

        mag_P_lat_in_DB = 20.0*np.log10(mag_P_lat_in)
        mag_CP_lat_in_DB = 20.0*np.log10(mag_CP_lat_in)

    #mag_CP_lat_out, phase, omegas = bode(P_lat_out*C_lat_out, plot=False, omega=omegas)


    print('Close window to end program')
    plt.show()

