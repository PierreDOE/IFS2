# -*- coding: utf-8 -*-
"""

Pierre DOERFLER & Batiste RIVIERE, January 2025

"""
from Flutter.flutter import Flutter
from Flutter.data import Data
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

if __name__ == "__main__":
    D = Data()
    F = Flutter(D.params)
    F.solve()
