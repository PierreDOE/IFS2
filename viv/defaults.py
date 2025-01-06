# -*- coding: utf-8 -*-
"""
Project name: galloping 
@author: Christophe Airiau
@data: novembre 2021

Inspired from Violette et al, JFS, 2012.

Defaults parameters

"""
import numpy as np


def default_parameters():
    """
    default parameters of the class viv
    """
    par = {}
    par["A"] = 12                   # coupling parameter
    par["mass_number"] = 0.02       # Mass number
    par["mu"] = 1                   # mass parameter
    par["gamma"] = 0.5              # fluid added damping coefficient
    par["xi"] = 1.                  # tension cable parameter
    par["nu"] = 0                   # 0: no damping, 1 : full case
    par["u"] = 1                    # reduced flow velocity
    par["epsilon"] = 1              # ???
    par["k_ref"] = 1                # wave number
    par["k_lim"] = [0, 2, 101]
    par["option"] = dict(verbose=True, plot=True)
    return par


