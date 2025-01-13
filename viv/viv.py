# -*- coding: utf-8 -*-
"""

 VIV analysis following Violette et al, 2012

 @author : christophe Airiau
 @date : october 2021


class VIV

"""
import numpy as np
import matplotlib.pyplot as plt
from viv.colored_messages import set_section

class VIV(object):
    """
    """
    def __init__(self, par):
        """
        initialisation
        """
        self.A = par["A"]                       # coupling parameter
        self.M = par["mass_number"]   # Mass number
        self.mu = par["mu"]                     # mass parameter
        self.gamma = par["gamma"]               # fluid added damping coefficient
        self.xi = par["xi"]                     # tension cable parameter
        self.nu = par["nu"]                     # 0: no damping, 1 : full case
        self.u = par["u"]                       # reduced flow velocity
        self.epsilon = par["epsilon"]
        self.option = par["option"]             # various options
        self.k_ref = par["k_ref"]
        self.k_lim = par["k_lim"]

        self.kinf = self.u * np.sqrt(1 - self.A * self.M)
        self.ksup = self.u * np.sqrt(1 + self.A * self.M)

        self.omega = None
        self.alpha = None
        self.delta = None
        self.k_table = None
        self.k = self.k_ref
        self.gain, self.phase = [], []
        self.coefs = None

    # ***********************************************
    #  SOLVER
    # ***********************************************
    def solve_omega(self):
        """
        return the roots of the dispersion equation
        """
        set_section("solve omega")
        self.set_k_table()
        self.omega = []
        for k in self.k_table:
            self.omega.append(sorted(np.roots(self.set_polynomial_coefs(k, self.u)), reverse=True))
            # print(f'k = {k}\tdet = {self.discriminant(k)}\tomega = {sorted(np.roots(self.set_polynomial_coefs(k, self.u)), reverse=True)}')


    def solve_gain(self):
        """
        G(omega) = \hat q / \hat y for omega in an interval
        return gain and phase list for all omega complex variables
        """
        for omg, k in zip(np.array(self.omega), self.k_table):
            gain = self.A * omg ** 2 / (omg ** 2 - self.u ** 2)
            self.gain.append(gain)
            self.phase.append(np.angle(-gain, deg=False))
            # print(f'omega = {omg} \t gain = {gain)} \t phase = {}')

    # ***********************************************
    #  GETTERS
    # ***********************************************
    def get_omega_table(self):
        """
        display eigen values on a nice table
        """
        test = self.u * np.sqrt(1 - (self.A * self.M) / 4) + 1j * self.u * np.sqrt(self.A * self.M) / 2
        print('valeurs omega_max à trouvé :', test)
        form = (" %+2.3f +i %+2.3f   " * 4)
        for omg, k in zip(np.array(self.omega[1:]), self.k_table):
            print(form % (omg[0].real, omg[0].imag, omg[1].real, omg[1].imag,
                          omg[2].real, omg[2].imag, omg[3].real, omg[3].imag))
            u = 0
            
            for i in omg:
                if np.abs(i - test) <= (1+1j) * self.epsilon 
            # r1 = np.sqrt(omg[0])
            # r2 = -np.sqrt(omg[0])
            # r3 = np.sqrt(omg[1])
            # r4 = -np.sqrt(omg[1])
            # print(form % (r1.real, r1.imag, r2.real, r2.imag,
            #               r3.real, r3.imag, r4.real, r4.imag))

    def get_parameters():
        """
        display parameters
        """
        print("alpha                   : %10.5f" % alpha)
        print("delta                   : %10.5f" % delta)
        print("beta                    : %10.5f" % beta)
        print("u                       : %10.5f" % self.u)

    # ***********************************************
    #  SETTERS
    # ***********************************************

    def set_damping_percentage(self, nu):
        """
        quantify the damping,  0 <= nu <= 1
        """
        if 0 <= nu <= 1:
            self.nu = nu
        else:
            self.nu = 0
        print("damping in %   : ", self.nu * 100)

    def set_k_table(self):
        """
        return the array of the non dimensional wave number k
        """
        self.k_table = np.linspace(self.k_lim[0], self.k_lim[1], self.k_lim[2])

    def set_alpha_beta_delta(self):
        """
        alpha  = A x M, beta = gamma / mu, delta = 1 - A x M
        """
        pass

    def determinant(self, omega, k):
        """
        Dispersion equation as a function of k and omega
        """
        pass

    def set_polynomial_coefs(self, k, u):
        """
        return the coefficient of the dispersion equation (determinant) 
        also found in self.coefs
        """
        self.coefs = [1,
                      0,
                      - k ** 2 - u ** 2 * (1 - self.A * self.M),
                      0,
                      k ** 2 * u ** 2]
        return self.coefs

    def discriminant(self,k):
        """
        discriminant with nu = 0
        """
        self.set_k_table()
        self.set_polynomial_coefs(k, self.u)
        return self.coefs[1] ** 2 - 4 * self.coefs[0] * self.coefs[2]

    def set_transfer_function(self, omega):
        """
        G(omega) = | \hat q / \hat y |,  omega real  
        """
        pass
        return abs(gain), np.angle(gain, deg=True)

    def set_omega_table(self, min_val=0., max_val=2.):
        """
        """
        omega = np.linspace(min_val, max_val, self.k_lim[2])
        pass
        

    # ***********************************************
    #  PLOTS 
    # ***********************************************
 
    def plot_k_omega_real(self):
        """

        """
        plt.figure()
        plt.title(r"Real part of $\omega$")
        plt.plot(self.k_table, np.abs(np.real(self.omega)), 'o')
        plt.xlabel(r"$k$")
        plt.ylabel(r'$Re(\omega)$')
        plt.grid()
        plt.show()


    def plot_k_omega_imag(self):
        """
        Complex part of omega for k
        """
        # roots = np.zeros_like(self.k_table)
        # u = 0
        # for i in np.abs(np.imag(self.omega)):
        #     for j in i:
        #         if j != any.roots[u]:
        #             roots[u] = j
        #     u += 1
        plt.figure()
        plt.title(r"Complexe part of $\omega$")
        plt.plot(self.k_table, np.imag(self.omega), 'o')
        plt.xlabel(r"$k$")
        plt.ylabel(r'$Im(\omega)$')
        plt.grid()
        plt.show()

    def plot_k_gain_modulus(self):
        """

        """
        plt.figure()
        plt.title(r"Real part of $\omega$")
        plt.plot(self.k_table, np.abs(self.gain), 'o')
        plt.xlabel(r"$k$")
        plt.ylabel(r'$Re(\omega)$')
        plt.grid()
        plt.show()

    def plot_k_gain_phase(self):
        """

        """
        plt.figure()
        plt.title(r"Real part of $\omega$")
        plt.plot(self.k_table, self.phase, 'o')
        plt.xlabel(r"$k$")
        plt.ylabel(r'$Re(\omega)$')
        plt.grid()
        plt.show()

