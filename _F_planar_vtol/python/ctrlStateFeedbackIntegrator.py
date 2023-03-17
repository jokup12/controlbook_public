#full state feedback
import numpy as np
import control as cnt
import VTOLParam as P


class ctrlStateFeedback:
    def __init__(self):
        trH = 3
        zetaH = .8
        trTheta = .3
        zetaTheta = .707
        trZ = 10.0*trTheta
        zetaZ = .707
        integrator_pole_lon = .1
        integrator_pole_lat = .1

        wnH = np.pi / (2 * trH * np.sqrt(1 - zetaH**2))
        #des_char_poly_H = (1, 2 * zetaH * wnH, wnH ** 2)
        #des_poles_lon = np.roots(des_char_poly_H)

        wnTheta = np.pi / (2 * trTheta * np.sqrt(1 - zetaTheta**2))
        #des_char_poly_theta = (1, 2 * zetaTheta * wnTheta, wnTheta ** 2)
        #des_poles_theta = np.roots(des_char_poly_theta)

        wnZ = np.pi / (2 * trZ * np.sqrt(1 - zetaZ**2))
        #des_char_poly_Z = (1, 2 * zetaZ * wnZ, wnZ ** 2)
        #des_poles_Z = np.roots(des_char_poly_Z)

        des_char_poly_lon = np.convolve([1, 2 * zetaH * wnH, wnH ** 2],
                                        [1, -integrator_pole_lon])
        des_poles_lon = np.roots(des_char_poly_lon)


        des_char_poly_lat = np.convolve(
                np.convolve([1, 2 * zetaZ * wnZ, wnZ**2],
                            [1, 2 * zetaTheta * wnTheta, wnTheta**2]),
                [1, -integrator_pole_lat])
        des_poles_lat = np.roots(des_char_poly_lat)


        #State Space Matrices
        Mass = P.mc + 2*P.mr

        A_lon = np.array([[0, 1], [0, 0]])
        B_lon = np.array([[1], [1/Mass]])
        C_lon = np.array([1, 0])
        D_lon = np.array([0])


        A1_lon = np.vstack((np.hstack((A_lon, np.zeros((np.size(A_lon,1),1)))),
                            np.hstack((-C_lon, np.array([0.0]))) ))
        print('A1_lon: ', A1_lon)
        B1_lon = np.vstack( (B_lon, 0.0) )
        print('B1_lon: ', B1_lon)



        A_lat = np.array([[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0],
                          [0.0, -P.Fe/Mass, -P.mu/Mass, 0.0], [0.0, 0.0, 0.0, 0.0]])
        B_lat = np.array([[0.0], [0.0], [0.0], [1/(P.Jc+2*P.mr*P.d**2)]])
        C_lat = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
        D_lat = np.array([[0.0], [0.0]])

        Cr_lat = C_lat[0]

        A1_lat = np.vstack((np.hstack((A_lat, np.zeros((np.size(A_lat, 1),1)))),
                            np.hstack((-Cr_lat, np.array([0.0]))) ))
        print('A1_lat: ', A1_lat)
        B1_lat = np.vstack((B_lat, 0.0))
        print('B1_lat: ', B1_lat)


        CC_lon = cnt.ctrb(A_lon, B_lon)
        CC_lat = cnt.ctrb(A_lat, B_lat)
        CC1_lon = cnt.ctrb(A1_lon, B1_lon)
        CC1_lat = cnt.ctrb(A1_lat, B1_lat)

        #print('rank LON: ', np.linalg.matrix_rank(CC_lon))
        #print('rank LAT: ', np.linalg.matrix_rank(CC_lat))

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(CC1_lon) != 3:
            print("The LON system is not controllable")
        else:
            K1_lon = cnt.place(A1_lon, B1_lon, des_poles_lon)
            self.K_lon = K1_lon[0][0:2]
            print('K_lon: ', self.K_lon)
            self.ki_lon = K1_lon[0][2]
            print('ki_lon: ', self.ki_lon)

        if np.linalg.matrix_rank(CC1_lat) != 5:
            print("The LAT system is not controllable")
        else:
            K1_lat = (cnt.place(A1_lat, B1_lat, des_poles_lat))
            self.K_lat = K1_lat[0][0:4]
            print('K_lat: ', self.K_lat)
            self.ki_lat = K1_lat[0][4]
            print('ki_lat: ', self.ki_lat)


        print('LON poles: ', des_poles_lon)
        print('LAT poles: ', des_poles_lat)

        self.h = 0.0
        self.hdot = 0.0
        self.hd1 = 0.0
        self.errorH = 0.0
        self.errorHd1 = 0.0
        self.integratorH = 0.0
        self.z = 0.0
        self.zdot = 0.0
        self.zd1 = 0.0
        self.errorZ = 0.0
        self.errorZd1 = 0.0
        self.integratorZ = 0.0
        self.theta = 0.0
        self.thetadot = 0.0
        self.thetad1 = 0.0
        self.errorTheta = 0.0
        self.errorThetad1 = 0.0
        self.integratorTheta = 0.0



    def update(self, ref, x):
        self.z = x[0][0]
        self.h = x[1][0]
        self.theta = x[2][0]
        #self.zdot = x[3][0]
        #self.hdot = x[4][0]
        #self.thetadot = x[5][0]



        z_ref = ref[0][0]
        h_ref = ref[1][0]

        self.errorH = h_ref - self.h
        self.errorZ = z_ref - self.z

        self.integratorH = self.integratorH \
            + (P.Ts / 2.0) * (self.errorH + self.errorHd1)
        self.integratorZ = self.integratorZ \
            + (P.Ts / 2.0) * (self.errorZ + self.errorZd1)

        #print('intH: ', self.integratorH)
        #print('intZ: ', self.integratorZ)


        self.zdot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts))*self.zdot \
                    + (2/(2*P.sigma + P.Ts))*(self.z - self.zd1)
        self.hdot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts))*self.hdot \
                    + (2/(2*P.sigma + P.Ts))*(self.h - self.hd1)
        self.thetadot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts))*self.thetadot \
                    + (2/(2*P.sigma + P.Ts))*(self.theta - self.thetad1)

        x_lon = np.array([[x[1, 0]], [x[4, 0]]])
        x_lat = np.array([ [x[0, 0]], [x[2, 0]], [x[3, 0]], [x[5, 0]] ])

        #F = self.kp * (self.ref - self.z) - self.kd * self.zdot
        F_tilde = -self.K_lon @ x_lon + self.ki_lon*self.integratorH
        F = P.Fe/np.cos(self.theta) + F_tilde[0]
        F = F.item(0)

        T = -self.K_lat @ x_lat + self.ki_lat*self.integratorZ
        T = T.item(0)

        print('F_t: ', F_tilde)
        #u = np.array([[F], [T]])
        #v = P.mixing @ np.array([[F], [T]])  # changes to fr and fl
        #self.fr = v[0][0]
        #self.fl = v[1][0]
        #self.fr = self.saturate(self.fr, P.Fmax)
        #self.fl = self.saturate(self.fl, P.Fmax)
        #v = np.array([[self.fr], [self.fl]])    #[fr, fl]T
        #u = np.array([[self.fr + self.fl], [P.d * (self.fr - self.fl)]])

        self.hd1 = self.h
        self.errorHd1 = self.errorH
        self.zd1 = self.z
        self.errorZd1 = self.errorZ
        self.thetad1 = self.theta
        self.errorThetad1 = self.errorTheta
        u = np.array([[F], [T]])
        return u

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u