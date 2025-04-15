'''
This is the "HelicalGear" class. It inherits from Component class.
This class is responsible for defining the properties of a HelicalGear object, such as:
I) Given Properties
--> 1) "name": a string of characters acting as a label
--> 2) "axis": a 3-element vector representing the axis along which the motor is rotating with respect to a defined reference frame
--> 3) "loc": a 3-element vector representing the location of the motor with respect to a defined reference frame
--> 4) "m_n": normal module of the gear expressed in [mm]
--> 5) "z": number of teeth of the gear
--> 6) "psi": helix angle expressed in [rad]
--> 7) "phi_n": normal pressure angle expressed in [rad]
--> 8) "Q_v": transmission accuracy
--> 9) "FW": tooth face width expressed in [mm]
--> 10) "material": a Material object representing the material properties of the gear
II) Calculated parameters
--> 1) "p_n": npmral circular pitch expressed in [mm]
--> 2) "p_t": transverse circular pitch expressed in [mm]
--> 3) "p_x": axial pitch expressed in [mm]
--> 4) "m_t": transverse module expressed in [mm]
--> 5) "d": pitch diameter expressed in [mm]
--> 6) "phi_t": transverse pressure angle expressed in [rad]
--> 7) "z_p": virtual number of teeth
--> 8) "phi_b": base helix angle expressed in [rad]
--> 9) "h_a": addendum height expressed in [mm]
--> 10) "h_f": dedendum height expressed in [mm]
--> 11) "h": tooth height expressed in [mm]
--> 12) "d_a": addendum diameter expressed in [mm]
--> 13) "d_f": dedendum diameter expressed in [mm]
III) Calculated loads
--> 1) "F_t": a 3-element vector representing the tangential force expressed in [N]
--> 2) "F_r": a 3-element vector representing the radial force expressed in [N]
--> 3) "F_a": a 3-element vector representing the axial force expressed in [N]
IV) Gear tooth verification parameters
--> 1) "sigma_max_fatigue": maximum gear tooth bending stress for fatigue expressed in [MPa]
--> 2) "power_source": a string of characters deining the power source type
--> 3) "K_0": overload factor
--> 4) "t_R": difference between dedendum radius and internal gear radius expressed in [mm]
--> 5) "m_B": ratio between t_R and tooth height
--> 6) "K_B": rim-thickness factor
--> 7) "K_v": dynamic factor
--> 8) "C_mc": coefficient 1 for load distribution factor
--> 9) "C_pf": coefficient 2 for load distribution factor
--> 10) "C_pm": coefficient 3 for load distribution factor
--> 11) "C_ma": coefficient 4 for load distribution factor
--> 12) "C_e":coefficient 5 for load distribution factor
--> 13) "Y": coefficient 6 for load distribution factor
--> 14) "K_H": load distribution factor
--> 15) "K_S": size factor
--> 16) "J_p": geometry factor
--> 17) "J_pp": modifying factor
--> 18) "Y_J": bending strength geometry factor
--> 19) "sigma_FP": bending fatigue strength expressed in [MPa]
--> 20) "Y_N": stress cycle life factor (bending)
--> 21) "Y_theta": temperature factor
--> 22) "Y_Z": reliability factor
--> 23) "bendingSF": bending safety factor
--> 24) "sigma_max_pitting": maximum gear contact (pitting resistance) stess expressed in [MPa]
--> 25) "Z_E": elastic coefficient expressed in [MPa^(1/2)]
--> 26) "Z_I": surface strength geometry
--> 27) "m_N": load sharing ratio
--> 28) "sigma_HP": contact strength expressed in [MPa]
--> 29) "Z_N": stress cycle life factor (wear)
--> 30) "Z_W": hardness-ratio factor
--> 31) "wearSF": wear safety factor
Moreover, it can carry out all the necessary calculations.
'''
import numpy as np
from .component import Component
from .force import Force
from math import pi, cos, sin, tan, atan, ceil, sqrt, log
from .makima2dInterpolator import Makima2DInterpolator
class HelicalGear(Component):

    # Overload factor reference matrix
    K0_ref = np.array([[1, 1.25, 1.75], [1.25, 1.5, 2], [1.5, 1.75, 2.25]])
    # Size factor reference graph
    z_ref = np.concatenate((np.arange(12, 23, 1), np.arange(24, 31, 2), np.array([34, 38, 43, 50, 60, 75, 100, 150, 300, 400])))
    Y_ref = np.array([0.245, 0.261, 0.277, 0.29, 0.296, 0.303, 0.309, 0.314, 0.322, 0.328, 0.331, 0.337, 0.346, 0.353, 0.359, 0.371, 0.384, 0.397, 0.409, 0.422, 0.435, 0.447, 0.46, 0.472, 0.48])
    # Bending strength geometry factor reference graphs
    psi_Jp_ref = np.arange(5, 35, 1)
    z_Jp_ref = np.array([20, 30, 60, 150, 500])
    Jp_ref = np.array([[0.465, 0.475, 0.48, 0.487, 0.492, 0.495, 0.497, 0.5, 0.502, 0.505, 0.506, 0.507, 0.508, 0.507, 0.506, 0.505, 0.502, 0.5, 0.497, 0.494, 0.49, 0.486, 0.48, 0.476, 0.471, 0.465, 0.458, 0.452, 0.445, 0.439], 
                       [0.525, 0.53, 0.535, 0.54, 0.542, 0.547, 0.55, 0.552, 0.553, 0.554, 0.555, 0.556, 0.555, 0.554, 0.552, 0.551, 0.549, 0.545, 0.54, 0.537, 0.532, 0.527, 0.52, 0.515, 0.507, 0.5, 0.492, 0.484, 0.476, 0.47], 
                       [0.58, 0.585, 0.595, 0.6, 0.602, 0.605, 0.61, 0.612, 0.614, 0.615, 0.615, 0.615, 0.614, 0.611, 0.608, 0.605, 0.6, 0.595, 0.59, 0.582, 0.577, 0.571, 0.56, 0.555, 0.545, 0.537, 0.527, 0.517, 0.507, 0.495], 
                       [0.62, 0.63, 0.635, 0.64, 0.645, 0.647, 0.652, 0.655, 0.657, 0.657, 0.656, 0.655, 0.654, 0.65, 0.645, 0.64, 0.635, 0.63, 0.622, 0.617, 0.61, 0.6, 0.592, 0.582, 0.575, 0.56, 0.552, 0.541, 0.53, 0.518], 
                       [0.65, 0.655, 0.66, 0.665, 0.67, 0.675, 0.677, 0.68, 0.681, 0.681, 0.681, 0.681, 0.679, 0.676, 0.673, 0.667, 0.66, 0.655, 0.647, 0.64, 0.631, 0.62, 0.613, 0.602, 0.592, 0.58, 0.57, 0.557, 0.545, 0.532]])
    psi_Jpp_ref = np.arange(5, 36, 5)
    z_Jpp_ref = np.array([20, 30, 50, 75, 150, 500])
    Jpp_ref = np.array([[0.927, 0.929, 0.93, 0.932, 0.938, 0.943, 0.953],
                        [0.952, 0.954, 0.957, 0.959, 0.961, 0.965, 0.973],
                        [0.98, 0.981, 0.981, 0.982, 0.983, 0.987, 0.992],
                        np.ones(7),
                        [1.02, 1.02, 1.02, 1.019, 1.017, 1.013, 1.01],
                        [1.038, 1.037, 1.035, 1.032, 1.03, 1.026, 1.02]])
    rel_ref = np.array([.9999, .999, .99, .9, .5])
    YZ_ref = np.array([1.5, 1.25, 1, 0.85, 0.7])
    ZE_ref = np.array([[191, 181, 179, 174, 162, 158], 
                       [181, 174, 172, 168, 158, 154], 
                       [179, 172, 170, 166, 156, 152], 
                       [174, 168, 166, 163, 154, 149], 
                       [162, 158, 156, 154, 145, 141], 
                       [158, 154, 152, 149, 141, 137]])

    # Constructor
    def __init__(self, name, axis, loc, m_n, z, psi, phi_n, Q_v, FW, material):
        # Given parameters
        super().__init__(name=name, material=material, axis=axis, loc=loc, F_tot=None, T_tot=None, omega=None)
        self.m_n = m_n
        self.z = z
        self.psi = psi * pi / 180
        self.phi_n = phi_n * pi / 180
        self.Q_v = Q_v
        self.FW = FW
        # Calculated parameters
        self.p_n = self.m_n * pi
        self.p_t = self.p_n / cos(self.psi)
        self.p_x = self.p_t / tan(self.psi)
        self.m_t = self.m_n / cos(self.psi)
        self.d = self.m_t * self.z
        self.phi_t = atan(tan(self.phi_n) / cos(self.psi))
        self.z_p = ceil(self.z / (cos(self.psi) ** 3))
        self.phi_b = atan(tan(self.psi) * cos(self.phi_n))
        self.h_a = self.m_n
        self.h_f = 1.25 * self.m_n # assumption
        self.h = self.h_a + self.h_f
        self.d_a = self.d + 2 * self.h_a
        self.d_f = self.d - 2 * self.h_f
    
    # Calculate Forces
    def calculateForces(self, mesh, shaft):
        if self.name == mesh.drivingGear.name:
            self.F_t = np.cross(self.T_tot.torque, 2 / self.d * mesh.axis) * 1e3
            mesh.drivingGear.F_t = self.F_t
            magF_t = sqrt(np.sum(self.F_t * self.F_t))
            self.F_r = -magF_t * tan(self.phi_n) / cos(self.psi) * mesh.axis
            mesh.drivingGear.F_r = self.F_r
            self.F_a = np.sign(np.sum(shaft.axis)) * np.sign(self.psi) * np.abs(np.cross(mesh.axis, self.F_t * tan(np.abs(self.psi))))
            mesh.drivingGear.F_a = self.F_a
            Floc = self.d / 2 * np.abs(mesh.axis) + self.loc
            self.F_tot = Force(self.F_t + self.F_r + self.F_a, Floc)
            mesh.drivingGear.F_tot = self.F_tot
        else:
            self.F_t = -mesh.drivingGear.F_t
            self.F_r = -mesh.drivingGear.F_r
            self.F_a = -mesh.drivingGear.F_a
            Floc = -self.d / 2 * np.abs(mesh.axis) + self.loc
            self.F_tot = Force(-mesh.drivingGear.F_tot.force, Floc)
            mesh.drivenGear.F_tot = self.F_tot
    
    # Maximum tooth gear bending stress equation for fatigue
    def calculateSigmaMaxFatigue(self, mesh, powerSource, drivenMachine, dShaft, Ce, teethCond, lShaft, useCond):
        fw_min = min(mesh.drivingGear.FW, mesh.drivenGear.FW)
        # Overload factor K_0
        if powerSource == "Uniform":
            i = 1
        elif powerSource == "Light shock":
            i = 2
        elif powerSource == "Medium shock":
            i = 3
        else:
            raise ValueError("Wrong input for power source")
        if drivenMachine == "Uniform":
            j = 1
        elif drivenMachine == "Moderate shock":
            j = 2
        elif drivenMachine == "Heavy shock":
            j = 3
        else:
            raise ValueError("Wrong input for driven machine")
        self.K_0 = self.__class__.K0_ref[i-1][j-1]
        # Rim-thickness factor
        self.t_R = (self.d_f - dShaft) / 2
        self.m_B = self.t_R / self.h
        if self.m_B < 1.2:
            self.K_B = 1.6 * log(2.242 / self.m_B)
        else:
            self.K_B = 1
        # Dynamic factor
        v = np.sum(np.abs(self.omega)) * self.d * 1e-3
        B = 0.25 * (12 - self.Q_v) ** (2/3)
        A = 50 + 56 * (1 - B)
        self.K_v = ((A + sqrt(200 * v)) / A) ** B
        # Load distribution factor
        self.C_e = Ce
        if teethCond == "uncrowned teeth":
            self.C_mc = 1
        elif teethCond == "crowned teeth":
            self.C_mc = 0.8
        else:
            raise ValueError("Invalid input for teeth condition.")
        if fw_min / 25.4 <= 1:
            self.C_pf = fw_min / 10 / self.d - 0.025
        elif fw_min / 25.4 > 1 and fw_min / 25.4 <= 17:
            self.C_pf = fw_min / 10 / self.d - 0.0375 + 0.0125 * fw_min / 25.4
        else:
            self.C_pf = fw_min / 10 / self.d - 0.1109 + 0.0207 * fw_min / 25.4 - 0.000228 * (fw_min / 25.4) ** 2
        S1 = abs(lShaft / 2 - np.sum(np.abs(self.loc)))
        if S1 / lShaft < 0.175:
            self.C_pm = 1
        else:
            self.C_pm = 1.1
        if useCond == "Open gearing":
            A = 0.247
            B = .0167
            C = -0.765 * 1e-4
        elif useCond == "Commercial, enclosed units":
            A = 0.127
            B = 0.0158
            C = -0.93 * 1e-4
        elif useCond == "Precision, enclosed units":
            A = 0.0675
            B = 0.0128
            C = -0.926 * 1e-4
        elif useCond == "Extraprecision enclosed gear units":
            A = 0.0036
            B = 0.0102
            C = -0.822 * 1e-4
        else:
            raise ValueError("wrong input for use.")
        self.C_ma = A + B * fw_min / 25.4 + C * (fw_min / 25.4) ** 2
        self.K_H = 1 + self.C_mc * (self.C_pf * self.C_pm + self.C_ma * self.C_e)
        # Size factor
        self.Y = np.interp(self.z_p, self.z_ref, self.Y_ref)
        self.K_S = 0.843 * (fw_min * self.m_t * sqrt(self.Y)) ** 0.0535
        # Bending strength geometry factor
        interp_func1 = Makima2DInterpolator(self.psi_Jp_ref, self.z_Jp_ref, self.Jp_ref)
        self.J_p = interp_func1(abs(self.psi * 180 / pi), self.z)
        interp_func2 = Makima2DInterpolator(self.psi_Jpp_ref, self.z_Jpp_ref, self.Jpp_ref)
        self.J_pp = interp_func2(abs(self.psi * 180 / pi), self.z)
        self.Y_J = self.J_p * self.J_pp
        self.sigma_max_fatigue = np.sum(np.abs(self.F_t)) * self.K_0 * self.K_B * self.K_v * self.K_H * self.K_S / fw_min / self.m_t / self.Y_J
    
    # Calulcate bending safety factor
    def calculateBendingSF(self, sigma_FP, b_YN, e_YN, N, temp, rel):
        self.sigma_FP = sigma_FP
        self.Y_N = b_YN * N ** e_YN
        if temp <= 120:
            self.Y_theta = 1
        if rel in self.__class__.rel_ref:
            print("I'm here!")
            ind = np.where(self.__class__.rel_ref == rel)
            self.Y_Z = self.__class__.YZ_ref[ind]
        else:
            self.Y_Z = np.interp(rel, self.__class__.rel_ref, self.__class__.YZ_ref)
        self.bendingSF = self.sigma_FP * self.Y_N / self.sigma_max_fatigue / self.Y_theta / self.Y_Z
    
    # Calculate sigma max pitting
    def calculateSigmaMaxPitting(self, mesh, Z_R):
        # Elastic coefficient
        if mesh.drivingGear.material.name == "Steel":
            i = 1
        elif mesh.drivingGear.material.name == "Malleable iron":
            i = 2
        elif mesh.drivingGear.material.name == "Nodular iron":
            i = 3
        elif mesh.drivingGear.material.name == "Cast iron":
            i = 4
        elif mesh.drivingGear.material.name == "Aluminum bronze":
            i = 5
        elif mesh.drivingGear.material.name == "Tin bronze":
            i = 6
        else:
            raise ValueError("Wrong input for pinion material.")
        if mesh.drivenGear.material.name == "Steel":
            j = 1
        elif mesh.drivenGear.material.name == "Malleable iron":
            j = 2
        elif mesh.drivenGear.material.name == "Nodular iron":
            j = 3
        elif mesh.drivenGear.material.name == "Cast iron":
            j = 4
        elif mesh.drivenGear.material.name == "Aluminum bronze":
            j = 5
        elif mesh.drivenGear.material.name == "Tin bronze":
            j = 6
        else:
            raise ValueError("Wrong input for gear material.")
        self.Z_E = self.__class__.ZE_ref[i-1][j-1]
        # Surface strength geometry
        pG = mesh.drivingGear
        gG = mesh.drivenGear
        f1_ZA = sqrt((pG.d / 2) ** 2 + pG.h_a ** 2)
        f2_Z = (pG.d + gG.d) * sin(self.phi_t) / 2
        Z_A = min([f1_ZA, f2_Z])
        f1_ZB = sqrt((gG.d / 2) ** 2 + gG.h_a ** 2)
        Z_B = min([f1_ZB, f2_Z])
        Z = Z_A + Z_B - f2_Z
        self.m_N = self.p_n * cos(self.phi_n) / 0.95 / Z
        if mesh.type == "External":
            self.Z_I = cos(self.phi_t) * sin(self.phi_t) * mesh.m_G / 2 / self.m_N / (mesh.m_G + 1)
        elif mesh.type == "Internal":
            self.Z_I = cos(self.phi_t) * sin(self.phi_t) * mesh.m_G / 2 / self.m_N / (mesh.m_G - 1)
        else:
            raise ValueError("Gear mesh type invalid.")
        self.sigma_max_pitting = self.Z_E * np.sqrt(np.sum(np.abs(self.F_t)) * self.K_0 * self.K_v * self.K_S * self.K_H * Z_R / self.FW / self.d / self.Z_I)
    # Calculate wear safety factor
    def calculateWearSF(self, sigma_HP, b_ZN, e_ZN, N, mesh):
        self.sigma_HP = sigma_HP
        self.Z_N = b_ZN * N ** e_ZN
        HB_ratio = mesh.drivingGear.material.HB / mesh.drivenGear.material.HB
        if HB_ratio < 1.2:
            A_p = 0
        elif HB_ratio >= 1.2 and HB_ratio <= 1.7:
            A_p = (8.98 * HB_ratio - 8.29) * 1e-3
        else:
            A_p = 0.00698
        self.Z_W = 1 + A_p * (mesh.m_G - 1)
        self.wearSF = self.sigma_HP * self.Z_N * self.Z_W / self.sigma_max_pitting / self.Y_theta / self.Y_Z