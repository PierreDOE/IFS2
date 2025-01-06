# -*- coding: utf-8 -*-
"""

 This program reproduces figs 5 --> 11 and 15 --> 17 from NACA 496
 Figs 12 --> 14 require multiple executions of dataset number 2 with
 appropriately modified input quantities
 Implemented using NACA 496 equations and 2DOF solution method
 Boyd Perry, III
 NASA-Langley Research Center
 May 2015
 class DataFigure
"""

from viv.colored_messages import *
import pandas as pd

class DataFigure(object):
    """
    """
    def __init__(self, par):
        """
        initialization
        """
        self.dataset = par["dataset"]
        self.case = par["case"]
        self.fig = par["fig"]
        self.range = par["range"]
        self.figtype = par["figtype"]
        self.title = par["title"]
        self.label = None
        self.var = None


def fill_data():
    """
    Three formats are returned, it is too much !

    it allows to add many more figures 
    Return:
         a list of dictionary
         a list of class DataFigure
         a pandas Frame
    """
    set_msg("fill the figure data")
    d = []
    par = []

    par.append(dict(index=0, fig=5, case=3, dataset=1, figtype="A", range=[0, 1, 0, 160],
                    title="Standard Case"))
    par.append(dict(index=1, fig=7, case=2, dataset=1, figtype="A", range=[0, 4, -0.008, 0.014],
                    title=r"Standard Case $x_\beta =$ %f"))
    par.append(dict(index=2, fig=9, case=1, dataset=1, figtype="A", range=[0, 10, -30, 120],
                    title="Standard Case"))
    par.append(dict(index=3, fig=6, case=3, dataset=1, figtype="B", range=[0, 180, 0, 22],
                    title="Standard Case"))
    par.append(dict(index=4, fig=8, case=2, dataset=1, figtype="B", range=[0, 0.014, 0, 1.5],
                    title=r"Standard Case $x_\beta = $  %f"))
    par.append(dict(index=5, fig=10, case=1, dataset=1, figtype="B", range=[0, 22, 0, 1.5],
                    title="Standard Case"))
    par.append(dict(index=6, fig=11, case=1, dataset=2, figtype="C", range=[0, 5./3., 0, 1.75],
                    title="Standard Case"))
    par.append(dict(index=7, fig=15, case=1, dataset=3, figtype="D", range=[0, 1.5, 0, 50],
                    title=r"Wing A for $x_\alpha =$ %f"))
    par.append(dict(index=8, fig=16, case=2, dataset=4, figtype="D", range=[0, 1.8, 0, 50],
                    title=r"Wing B for $x_\beta =$ %f"))
    par.append(dict(index=9, fig=17, case=3, dataset=4, figtype="D", range=[0, 2.8, 0, 40],
                    title=r"Wing B"))
    for dic in par:
        d.append(DataFigure(dic))
    return par, d, pd.DataFrame(par)
