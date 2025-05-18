'''
This is the GeometricStressRaiser class.
It contains the data and the necessary calculations via interpolation to get the final factor.
The properties are:
1) "Kt_B": bending geometric stress raiser
2) "Kt_N": normal geometric stress raiser
3) "Kt_T": torsion geometric stress raiser
'''
import numpy as np
from .makima2dInterpolator import Makima2DInterpolator
import matplotlib.pyplot as plt
class GeometricStressRaiser():

    # Bending: notch radius to smaller diameter ratio
    r2d_B = np.arange(0.01, 0.31, 0.01)
    # Bending: bigger diamaeter to smaller diamaeter ratio
    D2d_B = np.array([1.01, 1.1, 1.5, 3, 6])
    # Reference bending geometric stress raiser
    Kt_Bref = np.array([[2.1, 1.9, 1.7, 1.62, 1.56, 1.5, 1.45, 1.42, 1.39, 1.36, 1.33, 1.31, 1.29, 
                        1.27, 1.25, 1.23, 1.22, 1.21, 1.2, 1.19, 1.18, 1.17, 1.16, 1.16, 1.15, 
                        1.15, 1.15, 1.14, 1.14, 1.14], 
                        [2.7, 2.4, 2.18, 2.01, 1.9, 1.79, 1.72, 1.67, 1.63, 1.59, 1.55, 1.53, 1.5, 
                        1.47, 1.45, 1.43, 1.4, 1.39, 1.37, 1.35, 1.34, 1.32, 1.31, 1.3, 1.29, 
                        1.28, 1.27, 1.27, 1.27, 1.26],
                        [3, 2.7, 2.4, 2.22, 2.07, 1.95, 1.85, 1.79, 1.72, 1.67, 1.63, 1.6, 1.57, 
                        1.54, 1.5, 1.47, 1.45, 1.43, 1.41, 1.4, 1.38, 1.36, 1.35, 1.34, 1.33, 
                        1.31, 1.31, 1.3, 1.29, 1.29], 
                        [3.3, 3, 2.7, 2.44, 2.3, 2.17, 2.06, 1.98, 1.9, 1.83, 1.77, 1.71, 1.67, 
                        1.63, 1.57, 1.55, 1.53, 1.5, 1.48, 1.46, 1.44, 1.41, 1.4, 1.39, 1.36, 
                        1.35, 1.34, 1.33, 1.32, 1.31],
                        [3.6, 3.2, 2.8, 2.58, 2.4, 2.26, 2.15, 2.04, 1.96, 1.88, 1.83, 1.78, 1.72, 
                        1.68, 1.63, 1.6, 1.57, 1.54, 1.53, 1.5, 1.48, 1.45, 1.43, 1.41, 1.39, 
                        1.37, 1.35, 1.35, 1.34, 1.34]])
    # Bending: Makima Interpolation Function
    interp_func_B = Makima2DInterpolator(r2d_B, D2d_B, Kt_Bref)

    # Normal: notch radius to smaller diameter ratio
    r2d_N = np.arange(0.05, 0.31, 0.01)
    # Normal: bigger diameter to smaller diameter ratio
    D2d_N = np.array([1.01, 1.05, 1.2, 1.5, 2])
    # Reference normal geometric stress raiser
    Kt_Nref = np.array([[1.37, 1.33, 1.3, 1.27, 1.24, 1.22, 1.2, 1.19, 1.18, 1.17, 1.17,
                        1.16, 1.15, 1.15, 1.15, 1.14, 1.14, 1.14, 1.14, 1.14, 1.14, 1.14, 1.14,
                        1.14, 1.14, 1.14],
                        [1.74, 1.66, 1.6, 1.54, 1.5, 1.45, 1.42, 1.4, 1.39, 1.37, 1.36,
                        1.34, 1.33, 1.31, 1.3, 1.29, 1.29, 1.28, 1.28, 1.27, 1.26, 1.26, 1.25,
                        1.25, 1.25, 1.25],
                        [2.18, 2.02, 1.92, 1.84, 1.75, 1.7, 1.63, 1.6, 1.57, 1.54, 1.52,
                        1.5, 1.47, 1.46, 1.44, 1.42, 1.41, 1.4, 1.39, 1.37, 1.37, 1.36, 1.35, 
                        1.34, 1.33, 1.33],
                        [2.43, 2.24, 2.13, 2.04, 1.94, 1.89, 1.83, 1.79, 1.74, 1.71, 1.69,
                        1.65, 1.62, 1.6, 1.57, 1.55, 1.53, 1.52, 1.5, 1.49, 1.47, 1.46, 1.45,
                        1.44, 1.43, 1.43],
                        [2.62, 2.43, 2.28, 2.16, 2.08, 2, 1.93, 1.88, 1.83, 1.8, 1.77,
                        1.73, 1.7, 1.67, 1.64, 1.62, 1.61, 1.59, 1.57, 1.55, 1.54, 1.52, 1.51,
                        1.5, 1.49, 1.48]])
    # Normal: Makima Interpolation Function
    interp_func_N = Makima2DInterpolator(r2d_N, D2d_N, Kt_Nref)

    # Torsion: notch radius to smaller diameter ratio
    r2d_T = np.arange(0.01, 0.31, 0.01)
    # Torsion: bigger diameter to smalleter diameter ratio
    D2d_T = np.array([1.1, 1.2, 1.3])
    # Refernce torsion geometric stress raiser
    Kt_Tref = np.array([[1.84, 1.55, 1.41, 1.34, 1.28, 1.22, 1.19, 1.17, 1.16, 1.14, 
                        1.14, 1.13, 1.12, 1.11, 1.11, 1.11, 1.11, 1.11, 1.1, 1.1, 1.1, 1.1, 
                        1.09, 1.09, 1.09, 1.08, 1.08, 1.08, 1.08, 1.08], 
                        [2.38, 2, 1.81, 1.69, 1.6, 1.51, 1.44, 1.39, 1.36, 1.33, 
                        1.3, 1.28, 1.26, 1.24, 1.21, 1.21, 1.2, 1.19, 1.18, 1.17, 1.17, 1.17, 
                        1.15, 1.14, 1.14, 1.13, 1.13, 1.13, 1.13, 1.13], 
                        [2.69, 2.3, 2.01, 1.86, 1.74, 1.65, 1.58, 1.51, 1.48, 1.45, 
                        1.41, 1.39, 1.37, 1.35, 1.33, 1.31, 1.3, 1.29, 1.28, 1.26, 1.25, 1.24, 
                        1.23, 1.22, 1.21, 1.21, 1.21, 1.21, 1.21, 1.2]])
    # Torsion: Makima Interpolation Functin
    interp_func_T = Makima2DInterpolator(r2d_T, D2d_T, Kt_Tref)

    # Constructor
    def __init__(self, r2d=0, D2d=0):
        # Bending
        self.Kt_Breq = self.__class__.interp_func_B(r2d, D2d)
        if self.Kt_Breq < 1:
            self.Kt_Breq = 1
        # Normal
        self.Kt_Nreq = self.__class__.interp_func_N(r2d, D2d)
        if self.Kt_Nreq < 1:
            self.Kt_Nreq = 1
        # Torsion
        self.Kt_Treq = self.__class__.interp_func_T(r2d, D2d)
        if self.Kt_Treq < 1:
            self.Kt_Treq = 1
    
    # Plot Reference Geometric Stress Raiser
    def plotKtRef(self, r2d=0, D2d=0, Kt=0, title=""):
        plt.figure()
        for i in range(Kt.shape[0]):
            plt.plot(r2d, Kt[i], linewidth=1.5,label="D/d:" + str(D2d[i]))
        plt.xlabel("r/d [-]")
        plt.ylabel("Kt [-]")
        plt.title(title + " Geometric Stress Raiser")
        plt.legend()
        plt.grid()
        plt.show()
    
    # Plot Reference Bending Geometric Stress Raiser
    def plotKtBRef(self):
        self.plotKtRef(self.__class__.r2d_B, self.__class__.D2d_B, self.__class__.Kt_Bref, "Bending")
    
    # Plot Reference Normal Geometric Stress Raiser
    def plotKtNRef(self):
        self.plotKtRef(self.__class__.r2d_N, self.__class__.D2d_N, self.__class__.Kt_Nref, "Normal")
    
    # Plot Reference Torsion Geometric Stress Raiser
    def plotKtTRef(self):
        self.plotKtRef(self.__class__.r2d_T, self.__class__.D2d_T, self.__class__.Kt_Tref, "Torsion")