#!/usr/bin/env python3
import numpy as np
import array as arr
import csv
import matplotlib.pyplot as plt
import pdb
from uncertainties import ufloat
import csv
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from functions import *

def main():
    y = np.array([16814.6158214537,	427.130169405182,	161.604697431015,	104.487160095613,	89.2648800814352,	72.5740172397344])
    x = np.array([1002., 42457., 118850., 197600., 239300., 297000.])
    dy = np.array([0.96108157124974,	0.374165056907366,	0.849876543151193,	0.280440350667696,	0.328851196169095,	0.137606149534729])
    B, dB=linear_regression_B(x,1/y,dy)
    print('Linear regression only B: volt =',B,'+', dB)







if __name__ == '__main__':
    main()
