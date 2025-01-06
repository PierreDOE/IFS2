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
import matplotlib
matplotlib.use('TkAgg')

#  canevas proposé
# Ne décommenter que quand l'implémentation est effective

par = default_parameters()
v = VIV(par)
v.set_alpha_beta_delta()
v.set_damping_percentage(0.01)
v.solve_omega()
v.get_omega_table()
v.plot_k_omega_real()
# v.plot_k_omega_imag()

# v.solve_gain()
# v.plot_k_gain_modulus()
# v.plot_k_gain_phase()

# v.set_omega_table(min_val=0., max_val=2)
plt.show()
 