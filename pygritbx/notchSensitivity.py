'''
This is the NotchSensitivity class.
It contains the data and the necessary calculations via interpolation to get the final factor.
'''
import numpy as np
import matplotlib.pyplot as plt
from .makima2dInterpolator import Makima2DInterpolator
class NotchSensitivity:
    # Notch radius [mm]
    notchRadius = np.arange(0.25, 4.1, 0.25)
    # Ultimate Strength [MPa]
    sigma_u = np.array([.4, .7, 1, 1.4])*1e3
    # Notch Sensitivity
    qRef = np.array([[0.48, 0.57, 0.62, 0.65, 0.68, 0.7, 0.72, 0.73, 0.74, 0.755, 0.76, 0.77, 0.775, 0.78, 0.785, 0.79],
                    [0.595, 0.7, 0.75, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.845, 0.85, 0.855, 0.86, 0.865, 0.87],
                    [0.74, 0.81, 0.85, 0.87, 0.885, 0.892, 0.895, 0.897, 0.899, 0.903, 0.905, 0.908, 0.912, 0.915, 0.92, 0.922],
                    [0.85, 0.9, 0.92, 0.94, 0.946, 0.952, 0.96, 0.964, 0.97, 0.973, 0.975, 0.9755, 0.976, 0.9763, 0.9765, 0.9767]])
    # Makima Interpolation Function
    interp_func = Makima2DInterpolator(notchRadius, sigma_u, qRef)
    
    # Constructor
    def __init__(self, notchRadius=0, sigma_u=0):
        self.qReq = self.__class__.interp_func(notchRadius, sigma_u)
    
    # Plot Reference Notch Sensitivity
    def plotqRef(self):
        plt.figure
        for i in range(len(self.sigma_u)):
            plt.plot(self.notchRadius, self.qRef[i], linewidth=1.5, label="Su = " + str(self.sigma_u[i]) + "[MPa]")
        plt.xlabel("Notch Radius [mm]")
        plt.ylabel("q [-]")
        plt.title("Notch Sensitivity")
        plt.legend()
        plt.grid()
        plt.show()