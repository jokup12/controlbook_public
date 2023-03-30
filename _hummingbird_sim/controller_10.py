import numpy as np
import hummingbirdParam as P


class ctrlPID:
    def __init__(self):
        # tuning parameters for pitch
        tr_pitch = .7
        zeta_pitch = .707
        self.ki_pitch = .03
        # gain calculation
        b_theta = P.ellT / (P.m1 * P.ell1 ** 2 + P.m2 * P.ell2 ** 2 + P.J1y + P.J2y)
        # print('b_theta: ', b_theta)
        wn_pitch = np.pi / (2 * tr_pitch * np.sqrt(1 - zeta_pitch ** 2))
        a0_pitch = 0.0
        a1_pitch = 0.0
        alpha0_pitch = wn_pitch ** 2
        alpha1_pitch = 2.0 * zeta_pitch * wn_pitch
        self.kp_pitch = (alpha0_pitch - a0_pitch) / b_theta
        self.kd_pitch = (alpha1_pitch - a1_pitch) / b_theta
        # print gains to terminal
        print('kp_pitch: ', self.kp_pitch)
        print('ki_pitch: ', self.ki_pitch)
        print('kd_pitch: ', self.kd_pitch)
        # delayed variables
        self.theta_d1 = 0.0
        self.theta_dot = 0.0
        self.integrator_theta = 0.0
        self.error_theta_d1 = 0.0  # pitch error delayed by 1


        # tuning parameters for roll
        tr_roll = .3
        zeta_roll = .707
        self.ki_roll = 0.0
        # gain calculation
        b_phi = 1/P.J1x
        # print('b_phi: ', b_phi)
        wn_roll = np.pi / (2 * tr_roll * np.sqrt(1 - zeta_roll ** 2))
        a0_roll = 0.0
        a1_roll = 0.0
        alpha0_roll = wn_roll ** 2
        alpha1_roll = 2.0 * zeta_roll * wn_roll
        self.kp_roll = (alpha0_roll - a0_roll) / b_phi
        self.kd_roll = (alpha1_roll - a1_roll) / b_phi
        # print gains to terminal
        print('kp_roll: ', self.kp_roll)
        print('ki_roll: ', self.ki_roll)
        print('kd_roll: ', self.kd_roll)
        # delayed variables
        self.phi_d1 = 0.0
        self.phi_dot = 0.0
        self.integrator_phi = 0.0
        self.error_phi_d1 = 0.0  # roll error delayed by 1
        self.phi_ref = 0.0


        #tuning parameters for yaw/psi
        tr_yaw = 12*tr_roll
        zeta_yaw = .9
        self.ki_yaw = 0.0015
        # gain calculation
        b_psi = P.ellT/(P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J1y + P.J2y)
        #print('b_psi: ', b_psi)
        wn_yaw = np.pi / (2 * tr_yaw * np.sqrt(1 - zeta_yaw ** 2))
        a0_yaw = 0.0
        a1_yaw = 0.0
        alpha0_yaw = wn_yaw ** 2
        alpha1_yaw = 2.0 * zeta_yaw * wn_yaw
        self.kp_yaw = (alpha0_yaw - a0_yaw) / b_psi
        self.kd_yaw = (alpha1_yaw - a1_yaw) / b_psi
        # print gains to terminal
        print('kp_yaw: ', self.kp_yaw)
        print('ki_yaw: ', self.ki_yaw)
        print('kd_yaw: ', self.kd_yaw)
        # delayed variables
        self.psi_d1 = 0.0
        self.psi_dot = 0.0
        self.integrator_psi = 0.0
        self.error_psi_d1 = 0.0  # yaw error delayed by 1


    def update(self, r, y):
        phi_ref = r[0][0]
        theta_ref = r[1][0]
        psi_ref = r[2][0]
        phi = y[0][0]
        theta = y[1][0]
        psi = y[2][0]

        #theta/pitch calculations
        self.theta_dot = ((2 * P.sigma - P.Ts) / (2 * P.sigma + P.Ts)) * self.theta_dot \
                       + (2 / (2 * P.sigma + P.Ts)) * (theta - self.theta_d1)
        force_fl = self.kp_pitch * (theta_ref - theta) - self.kd_pitch*self.theta_dot + self.ki_pitch*self.integrator_theta   #pitch
        Fe = ((P.m1 * P.ell1 + P.m2 * P.ell2) * P.g / P.ellT) / (np.cos(phi) * np.cos(theta))
        error_theta = theta_ref - theta
        self.integrator_theta = self.integrator_theta + (P.Ts / 2) * (error_theta + self.error_theta_d1)
        force_unsat = force_fl + Fe
        force = saturate(force_unsat, -P.force_max, P.force_max)
        self.theta_d1 = theta
        self.error_theta_d1 = error_theta
        if self.ki_pitch != 0.0:
            self.integrator_theta = self.integrator_theta + P.Ts/self.ki_pitch*(force-force_unsat)

        #phi_ref calculation, outer loop
        self.psi_dot = ((2 * P.sigma - P.Ts) / (2 * P.sigma + P.Ts)) * self.psi_dot \
                      + (2 / (2 * P.sigma + P.Ts)) * (psi - self.psi_d1)
        phi_ref = self.kp_yaw*(psi_ref - psi) - self.kd_yaw*self.psi_dot + self.ki_yaw*self.integrator_psi  #outer loop, yaw
        error_psi = psi_ref - psi
        self.integrator_psi = self.integrator_psi + (P.Ts / 2) * (error_psi + self.error_psi_d1)
        self.psi_d1 = psi
        self.error_psi_d1 = error_psi



        #tau calculations, inner loop
        self.phi_dot = ((2 * P.sigma - P.Ts) / (2 * P.sigma + P.Ts)) * self.phi_dot \
                       + (2 / (2 * P.sigma + P.Ts)) * (phi - self.phi_d1)
        T_fl = self.kp_roll * (phi_ref - phi) - self.kd_roll*self.phi_dot #+ self.ki_roll*self.integrator_phi #inner loop, roll
        Te = 0
        error_phi = phi_ref - phi
        self.phi_d1 = phi
        self.error_phi_d1 = error_phi
        #self.integrator_phi = self.integrator_phi + (P.Ts / 2) * (error_phi + self.error_phi_d1)
        torque_unsat = T_fl + Te
        # convert force and torque to pwm signals
        pwm = np.array([[force + torque_unsat / P.d],  # u_left
                        [force - torque_unsat / P.d]]) / (2 * P.km)  # r_right
        pwm = saturate(pwm, 0, .7)

        v = 2*P.km*pwm  #fl and fr values after saturation
        fl = v[0][0]
        fr = v[1][0]
        v = np.array([[fr+fl], [P.d*(fl-fr)]])
        torque = v.item(1)

        #antiwindup for outer loop
        if self.ki_yaw != 0.0:
            self.integrator_psi = self.integrator_psi + P.Ts/self.ki_yaw*(torque - torque_unsat)

        #debug
        print('psi integrator:, ', self.integrator_psi)
        #print('fl: ', fl, 'fr: ', fr, 'Total Force: ', fl+fr)
        #print('theta diff: ', theta_ref - theta)
        #print('phi: ', phi, 'phi_d1: ', self.phi_d1)
        #print('phiref: ', phi_ref - phi)
        #print('phidot:', self.phi_dot)

        # update all delayed variables




        # return pwm plus reference signals
        return pwm, np.array([[phi_ref], [theta_ref], [psi_ref]]), v



def saturate(u, low_limit, up_limit):
    if isinstance(u, float) is True:
        if u > up_limit:
            u = up_limit
        if u < low_limit:
            u = low_limit
    else:
        for i in range(0, u.shape[0]):
            if u[i][0] > up_limit:
                u[i][0] = up_limit
            if u[i][0] < low_limit:
                u[i][0] = low_limit
    return u




