# -*- coding: utf-8 -*-
"""

 VIV analysis following Violette et al, 2012

 @author : christophe Airiau
 @date : october 2021

"""
import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.scimath import sqrt as csqrt
from viv.defaults import default_parameters
from viv.viv import VIV
from viv.tools import test_values
#import matplotlib
#matplotlib.use('TkAgg')

cas = 2

######### Cas u=1
if cas == 1:
    par = default_parameters()
    v = VIV(par)
    v.set_alpha_beta_delta()
    v.set_damping_percentage(0.02)
    v.solve_omega()
    v.get_omega_table()
    v.plot_k_omega_real()
    v.plot_k_omega_imag()
    v.solve_gain()
    v.plot_k_gain_modulus()
    v.plot_k_gain_phase()
    plt.show()

######## Cas k=1
if cas == 2:
    par = default_parameters()
    v = VIV(par)
    v.set_alpha_beta_delta()
    v.set_damping_percentage(0.02)

    v.solve_omega_for_u_range(info=True, print_='all')  # Calcule omega pour chaque u
    v.get_omega_u_inf_sup()

    v.plot_u_omega_real()
    v.plot_u_omega_imag()

    v.solve_gain_u()
    v.plot_u_gain_modulus()
    v.plot_u_gain_phase()

    plt.show()

######### Cas u=k=1
if cas == 3:
    par = default_parameters()
    par["u_lim"] = [1, 1, 101]
    v = VIV(par)
    v.set_alpha_beta_delta()
    v.set_damping_percentage(0.02)

    v.solve_gain_u_k()
    v.plot_u_k_gain_modulus()
    v.plot_u_k_gain_phase()

    plt.show()
