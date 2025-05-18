'''
This is the FatigueLimitCorrectorFactors class.
It contains the data and the necessary calculations via interpolation to get the final factor.
The properties are:
1) "Cs": size factor
2) "Cl": load factor
3) "Cf": surface finish factor
'''
import numpy as np
from scipy.interpolate import RegularGridInterpolator
class FatigueLimitCorrectorFactors:
    
    # CONSTANTS
    # Graph reference diameter - [mm]
    d_ref = np.arange(10, 201, 10)
    # Reference size factor - [-]
    Cs_ref = np.array([1, 0.93, 0.875, 0.845, 0.82, 0.8, 0.785, 0.773, 0.76, 0.753, 
                       0.745, 0.74, 0.73, 0.725, 0.72, 0.715, 0.71, 0.705, 0.7, 0.695])
    # Reference load factor - [-]
    Cl_ref = 0.7
    # Reference ultimate tensile strength - [MPa]
    sigma_u_ref = np.arange(200, 1501, 100)
    # Reference roughness - [micrometers]
    Ra_ref = np.array([0.8, 1.6, 3.2, 6.3, 10, 40, 160])
    # Reference surface finish factor - [-]
    Cf_ref = np.array([np.ones(14), 
                       [1, 0.98, 0.97, 0.96, 0.96, 0.96] + [0.95]*8,
                       [1, 0.97, 0.96, 0.94, 0.93, 0.93, 0.92] + [0.91]*7,
                       [0.98, 0.95, 0.93, 0.92, 0.91, 0.9, 0.89, 0.88, 0.87] + [0.86]*5,
                       [0.97, 0.95, 0.92, 0.88, 0.87, 0.84, 0.82, 0.81, 0.78, 0.77, 0.76] + [0.75]*3,
                       [0.96, 0.92, 0.88, 0.84, 0.8, 0.76, 0.73, 0.7, 0.68, 0.66, 0.65, 0.63, 0.63, 0.62],
                       [0.94, 0.89, 0.85, 0.8, 0.76, 0.72, 0.67, 0.65, 0.61, 0.58, 0.55, 0.52, 0.5, 0.48]])
    # Interpolation Function
    interp_func = RegularGridInterpolator((Ra_ref, sigma_u_ref), Cf_ref)
    
    # Constructor
    def __init__(self, section=None):
        if section == None:
            raise ValueError("Section cannot be an empty argument.")
        
        # Cf
        self.Cf_req = self.__class__.interp_func([[section.Ra, section.material.sigma_u]])[0]
        # Cl
        if section.sigma_m_N != 0 or section.sigma_a_N != 0:
            self.Cl_req = self.__class__.Cl_ref
        else:
            self.Cl_req = 0
        # Cs
        if section.sigma_m_Mb != 0 or section.sigma_a_Mb != 0 or section.sigma_m_Mt != 0 or section.sigma_a_Mt != 0:
            self.Cs_req = np.interp(section.d, self.__class__.d_ref, self.__class__.Cs_ref)
        else:
            self.Cs_req = 1
        # C
        self.C = self.Cf_req * self.Cl_req * self.Cs_req