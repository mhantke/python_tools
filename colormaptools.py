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

Nal = 4.
al = [0.,0/Nal,1/Nal,2/Nal,3/Nal,4/Nal]
cdict_jet_lightbg2 = {'red':   [(al[0],1.0,1.0),
                                (al[1],0.0,0.0),
                                (al[2],0.0,0.0),
                                (al[3],1.0,1.0),
                                (al[4],1.0,1.0),
                                (al[5],0.5,0.5)],
                     
                      'green':  [(al[0],1.0,1.0),
                                 (al[1],0.0,0.0),
                                 (al[2],1.0,1.0),
                                 (al[3],1.0,1.0),
                                 (al[4],0.0,0.0),
                                 (al[5],0.0,0.0)],
                     
                      'blue':  [(al[0],1.0,1.0),
                                (al[1],1.0,1.0),
                                (al[2],1.0,1.0),
                                (al[3],0.0,0.0),
                                (al[4],0.0,0.0),
                                (al[5],0.0,0.0)]}

cmaps = {"jet_lightbg2":matplotlib.colors.LinearSegmentedColormap('jet_lightbg2', cdict_jet_lightbg2, 1024)}

def make_colorbar(filename,Nx,Ny,colormap=cm.jet,orientation="vertical"):
    X,Y = meshgrid(arange(Nx),arange(Ny))
    if orientation == "vertical":
        C = -Y
    else:
        C = X
    imsave(filename,C,cmap=colormap)

def complex_array_to_rgb(X, theme='dark', rmax=None):
    '''Takes an array of complex number and converts it to an array of [r, g, b],
    where phase gives hue and saturaton/value are given by the absolute value.
    Especially for use with imshow for complex plots.'''
    absmax = rmax or abs(X).max()
    Y = zeros(X.shape + (3,), dtype='float')
    Y[..., 0] = angle(X) / (2 * pi) % 1
    if theme == 'light':
        Y[..., 1] = np.clip(np.abs(X) / absmax, 0, 1)
        Y[..., 2] = 1
    elif theme == 'dark':
        Y[..., 1] = 1
        Y[..., 2] = np.clip(np.abs(X) / absmax, 0, 1)
    Y = matplotlib.colors.hsv_to_rgb(Y)
    return Y
