# -*- coding: utf-8 -*-
"""

 VIV analysis following Violette et al, 2012

 @author : christophe Airiau
 @date : october 2021


class VIV

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from viv.colored_messages import set_section

class VIV(object):
    """
    """
    def __init__(self, par):
        """
        initialisation
        """
        self.A = par["A"]                       # coupling parameter
        self.M = par["mass_number"]             # Mass number
        self.mu = par["mu"]                     # mass parameter
        self.gamma = par["gamma"]               # fluid added damping coefficient
        self.xi = par["xi"]                     # tension cable parameter
        self.nu = par["nu"]                     # 0: no damping, 1 : full case
        self.u = par["u"]                       # reduced flow velocity
        self.epsilon = par["epsilon"]
        self.beta = self.gamma / self.mu
        self.option = par["option"]             # various options
        self.k_ref = par["k_ref"]
        self.k_lim = par["k_lim"]
        self.u_lim = par["u_lim"]

        self.kinf = self.u * (1 - np.sqrt(self.A * self.M))
        self.ksup = self.u * (1 + np.sqrt(self.A * self.M))

        self.uinf = self.k_ref / (1 + np.sqrt(self.A * self.M))
        self.usup = self.k_ref / (1 - np.sqrt(self.A * self.M))

        self.omega_lim = par["omega_lim"]

        self.omega = None
        self.u_k_table = None
        self.alpha = None
        self.delta = None
        self.k_table = None
        self.u_table = None
        self.k = self.k_ref
        self.gain, self.phase = [], []
        self.coefs = None

    # ***********************************************
    #  SOLVER
    # ***********************************************
    def solve_omega(self, info=False, print_='else'):
        """
        return the roots of the dispersion equation
        """
        set_section("solve omega")
        self.set_k_table()
        #self.set_u_table()
        self.omega = []
        i = 0
        for k in self.k_table:
            self.omega.append(sorted(np.roots(self.set_polynomial_coefs(k,self.u)), reverse=True))
            if print_ != 'all':
                cond = i % 5 == 0
            else:
                cond = True
            if info and cond:
                print(f'k = {k:.2f}\tdiscriminant = {self.discriminant(k):.2f}\t\tomega = {sorted(np.roots(self.set_polynomial_coefs(k, self.u)), reverse=True)}')
            i += 1

    def solve_omega_for_u_range(self, info=False, print_='else'):
        """
        Calcule omega pour une plage de valeurs de u.
        """
        set_section("solve omega")
        #self.set_k_table()
        self.set_u_table()
        #self.set_u_table()
        self.omega = []
        i = 0
        for u in self.u_table:
            self.omega.append(sorted(np.roots(self.set_polynomial_coefs(self.k_ref,u)), reverse=True))
            if print_ != 'all':
                cond = i % 5 == 0
            else:
                cond = True
            if info and cond:
                print(f'k = {k:.2f}\tdiscriminant = {self.discriminant(k):.2f}\t\tomega = {sorted(np.roots(self.set_polynomial_coefs(k, self.u)), reverse=True)}')
            i += 1

    def solve_gain(self):
        """
        G(omega) = \hat q / \hat y for omega in an interval
        return gain and phase list for all omega complex variables
        """
        for omg, k in zip(np.array(self.omega), self.k_table):
            gain = self.A * omg ** 2 / (omg ** 2 - self.u ** 2)
            self.gain.append(gain)
            self.phase.append(np.angle(-gain, deg=False))
    
    def solve_gain_u(self):
        """
        G(omega) = \hat q / \hat y for omega in an interval
        return gain and phase list for all omega complex variables
        """
        for omg, u in zip(np.array(self.omega), self.u_table):
            gain = self.A * omg ** 2 / (omg ** 2 - self.u ** 2)
            self.gain.append(gain)
            self.phase.append(np.angle(-gain, deg=False)/np.pi)

    def solve_gain_u_k(self):
        """
        G(omega) = \hat q / \hat y for omega in an interval
        return gain and phase list for all omega complex variables
        """
        self.set_u_k_table()
        for omg in (np.array(self.u_k_table)):
            gain = self.A * omg ** 2 / (omg ** 2 - self.u_lim[0] ** 2 + 1j*self.nu*self.epsilon*omg*self.u_lim[0])
            self.gain.append(gain)
            self.phase.append(np.angle(-gain, deg=False)/np.pi)

    # ***********************************************
    #  GETTERS
    # ***********************************************
    def get_omega_table(self):
        """
        display eigen values on a nice table
        """
        set_section("Get omega")
        test = self.u * np.sqrt(1 - (self.A * self.M) / 4) + 1j * self.u * np.sqrt(self.A * self.M) / 2
        print('valeurs omega_max à trouvé :', test)
        form = (" %2.3f %+2.3f i\t\t" * 4)
        closest_distance = float('inf')
        closest_omega = None
        closest_k = None
        for omg, k in zip(np.array(self.omega[1:]), self.k_table):
            print(form % (omg[0].real, omg[0].imag, omg[1].real, omg[1].imag,
                          omg[2].real, omg[2].imag, omg[3].real, omg[3].imag))
            for omega_val in omg:
                distance = abs(test - omega_val)
                if distance < closest_distance:
                    closest_distance = distance
                    closest_omega = omega_val
                    closest_k = k

        print(f"L'omega le plus proche de test est : {closest_omega.real:+.3f} +i {closest_omega.imag:+.3f}")
        print(f"Il est associé à k = {closest_k}")

    def get_omega_inf_sup(self):
        """
        Retourne les 4 racines de omega correspondant à la valeur de k dans self.k_table la plus proche de k_inf.
        """
        set_section("Get omega inf and sup")
        k_inf_idx = np.argmin(np.abs(self.k_table - self.kinf))
        k_sup_idx = np.argmin(np.abs(self.k_table - self.ksup))
        omega_for_k_inf = self.omega[k_inf_idx]
        omega_for_k_sup = self.omega[k_sup_idx]
        print(f"Les 4 racines de omega_inf pour k_inf={self.k_table[k_inf_idx]}\u2243{self.kinf} sont :")
        print(f"Omega inf : \u00b1 {omega_for_k_inf[0].real:.3f} \u00b1 i {omega_for_k_inf[0].imag:.3f}\n")
        print(f"Les 4 racines de omega_sup pour k_sup={self.k_table[k_sup_idx]}\u2243{self.ksup} sont :")
        print(f"Omega inf : \u00b1 {omega_for_k_sup[0].real:.3f} \u00b1 i {omega_for_k_sup[0].imag:.3f}")

        return omega_for_k_inf

    def get_omega_u_inf_sup(self):
        """
        Retourne les 4 racines de omega correspondant à la valeur de u dans self.u_table la plus proche de u_inf.
        """
        set_section("Get omega inf and sup")
        u_inf_idx = np.argmin(np.abs(self.u_table - self.uinf))
        u_sup_idx = np.argmin(np.abs(self.u_table - self.usup))
        omega_for_u_inf = self.omega[u_inf_idx]
        omega_for_u_sup = self.omega[u_sup_idx]
        print(f"Les 4 racines de omega_inf pour u_inf={self.u_table[u_inf_idx]}\u2243{self.uinf} sont :")
        print(f"Omega inf : \u00b1 {omega_for_u_inf[0].real:.3f} \u00b1 i {omega_for_u_inf[0].imag:.3f}\n")
        print(f"Les 4 racines de omega_sup pour u_sup={self.u_table[u_sup_idx]}\u2243{self.usup} sont :")
        print(f"Omega inf : \u00b1 {omega_for_u_sup[0].real:.3f} \u00b1 i {omega_for_u_sup[0].imag:.3f}")

        return omega_for_u_inf

    def get_parameters(self):
        """
        display parameters
        """
        print("alpha                   : %10.5f" % self.alpha)
        print("delta                   : %10.5f" % self.delta)
        print("beta                    : %10.5f" % self.beta)
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
    
    def set_u_k_table(self):
        """
        return the array of the non dimensional wave number for u and k fixed
        """
        self.u_k_table = np.linspace(self.omega_lim[0], self.omega_lim[1], self.omega_lim[2])

    def set_u_table(self):
        """
        return the array of the non dimensional wave number k
        """
        self.u_table = np.linspace(self.u_lim[0], self.u_lim[1], self.u_lim[2])

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
        if self.nu == 0:
            self.coefs = [1,
                          0,
                          - k ** 2 - u ** 2 * (1 - self.A * self.M),
                          0,
                          k ** 2 * u ** 2]
        else:
            self.coefs = [1,
                          1j * self.nu * (-self.xi + (self.epsilon - self.beta) * u),
                          - k ** 2 - u ** 2 * (1 - self.A * self.M) + self.nu ** 2 * self.epsilon * (
                                      self.beta * u ** 2 + self.xi * u),
                          1j * self.nu * (self.beta * u ** 3 - self.epsilon * k ** 2 * self.u + self.xi * u ** 2),
                          k ** 2 * u ** 2]
        return self.coefs

    def discriminant(self,k):
        """
        discriminant with nu = 0
        """
        self.set_k_table()
        self.set_polynomial_coefs(k, self.u)
        return self.coefs[2] ** 2 - 4 * self.coefs[0] * self.coefs[4]

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
    def plot_u_omega_real(self):
        """Real part of omega for u
        """
        plt.figure()
        plt.title(f"Real part of $\\omega$ for $\\nu={self.nu}$\nRivière & Doerfler", usetex=False)
        plt.plot(self.u_table, np.abs(np.real(self.omega)), 'o')
        plt.xlabel(r"$u$")
        plt.ylabel(r'$Re(\omega)$')
        plt.grid()
        plt.savefig('Re_omega_u.png')
        plt.show()

    def plot_u_omega_imag(self):
        """
        Complex part of omega for u
        """
        plt.figure()
        plt.title(f"Complexe part of $\\omega$ for $\\nu={self.nu}$\nRivière & Doerfler", usetex=False)
        plt.plot(self.u_table, np.imag(self.omega), 'o')
        plt.xlabel(r"$u$")
        plt.ylabel(r'$Im(\omega)$')
        plt.grid()
        plt.savefig('Im_omega_u.png')
        plt.show()

    def plot_u_gain_modulus(self):
        """

        """
        plt.figure()
        plt.title(f"Modulus of $G$ for $\\nu={self.nu}$\nRivière & Doerfler", usetex=False)
        plt.plot(self.u_table, np.abs(self.gain), 'o')
        plt.ylim([0,80])
        plt.xlabel(r"$u$")
        plt.ylabel(r'$|G(\omega)|$')
        plt.grid()
        plt.savefig('modulus_omega_u')
        plt.show()

    def plot_u_gain_phase(self):
        """

        """
        plt.figure()
        plt.title(f"Phase of $G$ for $\\nu={self.nu}$\nRivière & Doerfler", usetex=False)
        plt.plot(self.u_table, self.phase, 'o')
        plt.xlabel(r"$u$")
        plt.ylabel(r'$\phi_G/\pi$')
        plt.grid()
        plt.savefig('Gain_omega_u')
        plt.show()

    def plot_k_omega_real(self):
        """

        """
        plt.figure()
        plt.title(f"Real part of $\\omega$ for $\\nu={self.nu}$\nRivière & Doerfler", usetex=False)
        plt.plot(self.k_table, np.abs(np.real(self.omega)), 'o')
        plt.xlabel(r"$k$")
        plt.ylabel(r'$Re(\omega)$')
        plt.grid()
        plt.savefig('Re_omega_k.png')
        plt.show()

    def plot_k_omega_imag(self):
        """
        Complex part of omega for k
        """
        plt.figure()
        plt.title(f"Complexe part of $\\omega$ for $\\nu={self.nu}$\nRivière & Doerfler", usetex=False)
        plt.plot(self.k_table, np.imag(self.omega), 'o')
        plt.xlabel(r"$k$")
        plt.ylabel(r'$Im(\omega)$')
        plt.grid()
        plt.savefig('Im_omega_k.png')
        plt.show()

    def plot_k_gain_modulus(self):
        """

        """
        plt.figure()
        plt.title(f"Modulus of $G$ for $\\nu={self.nu}$\nRivière & Doerfler", usetex=False)
        plt.plot(self.k_table, np.abs(self.gain), 'o')
        plt.xlabel(r"$k$")
        plt.ylabel(r'$|G(\omega)|$')
        plt.grid()
        plt.savefig('modulus_omega_k')
        plt.show()

    def plot_k_gain_phase(self):
        """

        """
        plt.figure()
        plt.title(f"Phase of $G$ for $\\nu={self.nu}$\nRivière & Doerfler", usetex=False)
        plt.plot(self.k_table, self.phase, 'o')
        plt.xlabel(r"$k$")
        plt.ylabel(r'$\phi_G/\pi$')
        plt.grid()
        plt.savefig('Gain_omega_k')
        plt.show()

    def plot_u_k_gain_modulus(self):
        """

        """
        plt.figure()
        plt.title(f"Modulus of $G$ for $\\nu={self.nu}$\nRivière & Doerfler", usetex=False)
        plt.plot(self.u_k_table, np.abs(self.gain), 'o')
        plt.xlabel(r"$omega$")
        plt.ylabel(r'$|G(\omega)|$')
        plt.grid()
        plt.savefig('modulus_omega_u_k')
        plt.show()

    def plot_u_k_gain_phase(self):
        """

        """
        plt.figure()
        plt.title(f"Phase of $G$ for $\\nu={self.nu}$\nRivière & Doerfler", usetex=False)
        plt.plot(self.u_k_table, self.phase, 'o')
        plt.xlabel(r"$omega$")
        plt.ylabel(r'$\phi_G/\pi$')
        plt.grid()
        plt.savefig('Gain_omega_u_k')
        plt.show()
