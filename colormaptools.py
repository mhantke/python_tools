#====================#
# Python tools - CXI #
#====================# 
#
# Author: Max Hantke
# Email: maxhantke@gmail.com


from matplotlib import colors
from pylab import *

cmaps = {}

cdict_jet = {'red':   [(0.0,0.0,0.0),
                       (0.375,0.0,0.0),
                       (0.625,1.0,1.0),
                       (0.875,1.0,1.0),
                       (1.0,0.5,0.5)],
             
             'green':  [(0.0,0.0,0.0),
                        (0.125,0.0,0.0),
                        (0.375,1.0,1.0),
                        (0.625,1.0,1.0),
                        (0.875,0.0,0.0),
                        (1.0,0.0,0.0)],
             
             'blue':  [(0.0,0.5,0.5),
                       (0.125,1.0,1.0),
                       (0.375,1.0,1.0),
                       (0.625,0.0,0.0),
                       (1.0,0.0,0.0)]}

cmaps = {"myjet":matplotlib.colors.LinearSegmentedColormap('myjet', cdict_jet, 1024)}

cdict_jet_lightbg = {'red':   [(0.0,1.0,1.0),
                               (0.125,0.0,0.0),
                               (0.375,0.0,0.0),
                               (0.625,1.0,1.0),
                               (0.875,1.0,1.0),
                               (1.0,0.5,0.5)],
                     
                     'green':  [(0.0,1.0,1.0),
                                (0.125,0.0,0.0),
                                (0.375,1.0,1.0),
                                (0.625,1.0,1.0),
                                (0.875,0.0,0.0),
                                (1.0,0.0,0.0)],
                     
                     'blue':  [(0.0,1.0,1.0),
                               (0.125,1.0,1.0),
                               (0.375,1.0,1.0),
                               (0.625,0.0,0.0),
                               (1.0,0.0,0.0)]}

cmaps = {"jet_lightbg":matplotlib.colors.LinearSegmentedColormap('jet_lightbg', cdict_jet_lightbg, 1024)}
