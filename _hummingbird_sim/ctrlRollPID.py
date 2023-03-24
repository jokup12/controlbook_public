import numpy as np
import hummingbirdParam as P


class ctrlRollPID:
    def __init__(self):
        # tuning parameters
        tr_roll = .3
        zeta_roll = .7
        self.ki_roll = 0
        # gain calculation
        b_phi = P.ellT/(P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J1y + P.J2y)
        #print('b_phi: ', b_phi)
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
        # sample rate of the controller
        self.Ts = P.Ts
        # dirty derivative parameters
        sigma = 0.05  # cutoff freq for dirty derivative
        self.beta = (2 * sigma - self.Ts) / (2 * sigma + self.Ts)
        # delayed variables
        self.phi_d1 = 0.
        self.phi_dot = 0.
        self.integrator_phi = 0.
        self.error_phi_d1 = 0.  # roll error delayed by 1

    def update(self, r, y):
        phi_ref = r[0][0]
        phi = y[1][0]
        force_fl = self.kp_roll*(phi_ref-phi) - self.kd_roll*self.phi_dot
        Fe = np.cos(phi)*(P.m1*P.ell1 + P.m2*P.ell2)*P.g/P.ellT
        # compute errors
        error_phi = phi_ref - phi
        # update differentiators
        self.phi_dot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts)) * self.phi_dot \
                + (2/(2*P.sigma + P.Ts))*(phi-self.phi_d1)
        
        # update integrators
        self.integrator_phi = self.integrator_phi + (P.Ts/2)*(error_phi + self.error_phi_d1)
        
        # roll control
        force_unsat = force_fl + Fe
        force = saturate(force_unsat, -P.force_max, P.force_max)
        torque = 0.
        # convert force and torque to pwm signals
        pwm = np.array([[force + torque / P.d],               # u_left
                      [force - torque / P.d]]) / (2 * P.km)   # r_right          
        pwm = saturate(pwm, 0, 1)
        # update all delayed variables
        self.phi_d1 = phi
        self.error_phi_d1 = error_phi
        # return pwm plus reference signals
        return pwm, np.array([[phi_ref], [0.], [0.]])


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




