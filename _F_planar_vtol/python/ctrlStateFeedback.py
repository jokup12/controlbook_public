#full state feedback
import numpy as np
import control as cnt
import VTOLParam as P


class ctrlStateFeedback:
    def __init__(self):
        trH = 3
        zetaH = .99
        trTheta = .2
        zetaTheta = .707
        trZ = 10.0*trTheta
        zetaZ = .707

        #wnH = np.pi / (2 * trH * np.sqrt(1 - zetaH**2))
        wnH = 2.2/trH
        des_char_poly_H = (1, 2 * zetaH * wnH, wnH ** 2)
        des_poles_lon = np.roots(des_char_poly_H)

        #wnTheta = np.pi / (2 * trTheta * np.sqrt(1 - zetaTheta**2))
        wnTheta = 2.2/trTheta
        des_char_poly_theta = (1, 2 * zetaTheta * wnTheta, wnTheta ** 2)
        des_poles_theta = np.roots(des_char_poly_theta)

        #wnZ = np.pi / (2 * trZ * np.sqrt(1 - zetaZ ** 2))
        wnZ = 2.2/trZ
        des_char_poly_Z = (1, 2 * zetaZ * wnZ, wnZ ** 2)
        des_poles_Z = np.roots(des_char_poly_Z)

        #des_poles_lat = np.array([des_poles_theta[0], des_poles_theta[1], des_poles_Z[0], des_poles_Z[1]])
        des_char_poly_lat = np.convolve([1, 2 * zetaZ * wnZ, wnZ ** 2],
                                         [1, 2 * zetaTheta * wnTheta, wnTheta ** 2])
        des_poles_lat = np.roots(des_char_poly_lat)

        #State Space Matrices
        Mass = P.mc + 2*P.mr

        A_lon = np.array([[0, 1], [0, 0]])
        B_lon = np.array([[1], [1/Mass]])
        C_lon = np.array([1, 0])
        D_lon = np.array([0])

        A_lat = np.array([[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0],
                          [0.0, -P.Fe/Mass, -P.mu/Mass, 0.0], [0.0, 0.0, 0.0, 0.0]])
        B_lat = np.array([[0.0], [0.0], [0.0], [1/(P.Jc+2*P.mr*P.d**2)]])
        C_lat = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
        D_lat = np.array([[0.0], [0.0]])


        CC_lon = cnt.ctrb(A_lon, B_lon)
        CC_lat = cnt.ctrb(A_lat, B_lat)
        print('rank LON: ', np.linalg.matrix_rank(CC_lon))
        print('rank LAT: ', np.linalg.matrix_rank(CC_lat))

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(CC_lon) != 2:
            print("The LON system is not controllable")
        else:
            self.K_lon = (cnt.place(A_lon, B_lon, des_poles_lon))
            Cr_lon = C_lon[0]
            self.kr_lon = -1.0 / (C_lon @ np.linalg.inv(A_lon - B_lon @ self.K_lon) @ B_lon)
        if np.linalg.matrix_rank(CC_lat) != 4:
            print("The LAT system is not controllable")
        else:
            self.K_lat = (cnt.place(A_lat, B_lat, des_poles_lat))
            Cr_lat = C_lat[0]
            self.kr_lat = -1.0 / (Cr_lat @ np.linalg.inv(A_lat - B_lat @ self.K_lat) @ B_lat)

        print('K_lon: ', self.K_lon)
        print('kr_lon: ', self.kr_lon)
        print('LON poles: ', des_poles_lon)
        print('K_lat: ', self.K_lat)
        print('kr_lat: ', self.kr_lat)
        print('LAT poles: ', des_poles_lat)

        self.h = 0.0
        self.hdot = 0.0
        self.hd1 = 0.0
        self.z = 0.0
        self.zdot = 0.0
        self.zd1 = 0.0
        self.theta = 0.0
        self.thetadot = 0.0
        self.thetad1 = 0.0



    def update(self, ref, x):
        self.z = x[0][0]
        self.h = x[1][0]
        self.theta = x[2][0]
        self.zdot = x[3][0]
        self.hdot = x[4][0]
        self.thetadot = x[5][0]

        z_ref = ref[0][0]
        h_ref = ref[1][0]



        #self.zdot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts))*self.zdot \
        #            + (2/(2*P.sigma + P.Ts))*(self.z - self.zd1)
        #self.hdot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts))*self.hdot \
        #            + (2/(2*P.sigma + P.Ts))*(self.h - self.hd1)
        #self.thetadot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts))*self.thetadot \
        #            + (2/(2*P.sigma + P.Ts))*(self.theta - self.thetad1)

        x_lon = np.array([[x[1, 0]], [x[4, 0]]])
        x_lat = np.array([ [x[0, 0]], [x[2, 0]], [x[3, 0]], [x[5, 0]] ])

        #F = self.kp * (self.ref - self.z) - self.kd * self.zdot
        F_tilde = -self.K_lon @ x_lon + self.kr_lon*h_ref
        F = P.Fe/np.cos(self.theta) + F_tilde[0,0]
        F = F.item(0)

        T = -self.K_lat @ x_lat + self.kr_lat*z_ref
        T = T.item(0)


        #u = np.array([[F], [T]])
        #v = P.mixing @ np.array([[F], [T]])  # changes to fr and fl
        #self.fr = v[0][0]
        #self.fl = v[1][0]
        #self.fr = self.saturate(self.fr, P.Fmax)
        #self.fl = self.saturate(self.fl, P.Fmax)
        #v = np.array([[self.fr], [self.fl]])    #[fr, fl]T
        #u = np.array([[self.fr + self.fl], [P.d * (self.fr - self.fl)]])

        self.hd1 = self.h
        self.zd1 = self.z
        self.thetad1 = self.theta
        u = np.array([[F], [T]])
        return u

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u