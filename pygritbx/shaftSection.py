'''
This is the "ShaftSection" class.
It's used to define certain sections of a "Shaft" component along its profile.
It has the following properties:
I) Given Properties:
--> 1) "name": a string of characters acting as a label
--> 2) "loc": a 3-D element vector representing the location of the section along the shaft
--> 3) "d": a scalar representing the diameter of the shaft at the specified cross-section expressed in [mm]
--> 4) "Ra": a scalar representing the surface roughness at the specified cross-section expressed in [micrometers]
--> 5) "material": a Material object representing the material of the shaft
II) Calculated Properties:
--> 1) "sigma_m_Mb": mean bending stress [MPa]
--> 2) "sigma_a_Mb": alternating bending stress [MPa]
--> 3) "sigma_m_N": mean normal stress [MPa]
--> 4) "sigma_a_N": alternating normal stress [MPa]
--> 5) "tau_m_Mt": mean torsional stress [MPa]
--> 6) "tau_a_Mt": alternating torsional stress [MPa]
--> 7) "Kt_B": bending stress concentration rasiser
--> 8) "Kt_N": normal stress concentration raiser
--> 9) "Kt_T": torsion stress concentration raiser
--> 10) "Kf_B": bending stress concentration factor
--> 11) "Kf_N": normal stress concentration factor
--> 12) "Kf_T": torsion stress concentration factor
'''
from .notchSensitivity import NotchSensitivity
from .fatigueLimitCorrectorFactors import FatigueLimitCorrectorFactors
from .geometricStressRaiser import GeometricStressRaiser
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
class ShaftSection:

    # Constructor
    def __init__(self, name="", loc=0, d=0, Ra=0, material=None):
        self.name = name
        self.loc = loc
        self.d = d
        self.Ra = Ra
        self.material = material
        # initialize mean and alternating stresses to 0
        self.sigma_m_Mb = 0
        self.sigma_a_Mb = 0
        self.sigma_m_N = 0
        self.sigma_a_N = 0
        self.tau_m_Mt = 0
        self.tau_a_Mt = 0
        # initialize stress concentration raisers
        self.Kt_B = 0
        self.Kt_N = 0
        self.Kt_T = 0
        # initialize stress concentration factor
        self.Kf_N = 1
        self.Kf_B = 1
        self.Kf_T = 1
        # initialize static safety factor
        self.staticSF = 0
    
    # Append fatigue stress intensification factor
    def appendKf(self, Kf=0, loadType=""):
        for kf, lt in zip(Kf, loadType):
            if lt == "Normal":
                self.Kf_N = self.Kf_N * kf
            elif lt == "Bending":
                self.Kf_B = self.Kf_B * kf
            elif lt == "Torsion":
                self.Kf_T = self.Kf_T * kf
            else:
                raise ValueError(f"'{lt}' is Invalid")
            
    # Calculate fatigue stress intensification factor
    def calculateFatigueIntensificationFactor(self):
        self.__class__.appendKf(self, [1 + self.q.qReq * (self.Kt_B - 1)], ["Bending"])
        self.__class__.appendKf(self, [1 + self.q.qReq * (self.Kt_N - 1)], ["Normal"])
        self.__class__.appendKf(self, [1 + self.q.qReq * (self.Kt_T - 1)], ["Torsion"])

    # Notch Sensitivity
    def addNotchSensitivity(self, notchRadius=0, sigma_u=0):
        self.q = NotchSensitivity(notchRadius=notchRadius, sigma_u=sigma_u)
        if self.Kt_B == None and self.Kt_N == None and self.Kt_T == None:
            self.__class__.calculateFatigueIntensificationFactor(self)
    
    #Geomtric Stress Raiser
    def addGeometricStressRaiser(self, r2d=0, D2d=0):
        temp = GeometricStressRaiser(r2d, D2d)
        self.Kt_B = temp.Kt_Breq
        self.Kt_N = temp.Kt_Nreq
        self.Kt_T = temp.Kt_Treq
        self.calculateFatigueIntensificationFactor()
    
    # Add fatigue limit corrector factors
    def addFLCF(self):
        self.FLCF = FatigueLimitCorrectorFactors(self)
        self.material.ComponentFatigueLimit(self)
    
    # Calculate Equivalent Stress
    def calculateSectionEquivalentStress(self):
        self.sigma_a_eq = sqrt((self.Kf_B * self.sigma_a_Mb + self.Kf_N * self.sigma_a_N/0.85) ** 2 + 3 * (self.Kf_T * self.tau_a_Mt) ** 2)
        self.sigma_m_eq = sqrt((self.Kf_B * self.sigma_m_Mb + self.Kf_N * self.sigma_m_N/0.85) ** 2 + 3 * (self.Kf_T * self.tau_m_Mt) ** 2)

    # Plot Haigh Diagram
    def plotHaighDiagram(self):
        # Calculate coordinates of specific points
        coeff1 = np.polyfit(np.array([0, self.material.sigma_y]), np.array([self.material.sigma_y, 0]), deg=1)
        coeff2 = np.polyfit(np.array([0, self.material.sigma_u]), np.array([self.material.sigma_Dm1C, 0]), deg=1)
        inter_x = (coeff2[1] - coeff1[1]) / (coeff1[0] - coeff2[0])
        inter_y = np.polyval(coeff1, inter_x)
        coeff3 = np.polyfit(np.array([0, self.sigma_m_eq]), np.array([0, self.sigma_a_eq]), deg=1)
        P_x = (coeff2[1] - coeff3[1]) / (coeff3[0] - coeff2[0])
        if P_x > inter_x:
            P_x = (coeff1[1] - coeff3[1]) / (coeff3[0] - coeff1[0])
        P_y = np.polyval(coeff3, P_x)
        # Plot
        figure, ax = plt.subplots()
        ax.plot(np.array([0, inter_x, self.material.sigma_y]), np.array([self.material.sigma_Dm1C, inter_y, 0]), 'k', linewidth=1.5)
        ax.plot(self.sigma_m_eq, self.sigma_a_eq, 'm*', linewidth=1.5)
        ax.plot(np.array([0, inter_x, self.material.sigma_u]), np.array([self.material.sigma_y, inter_y, 0]), 'k--', linewidth=1.5)
        ax.plot(0, self.material.sigma_Dm1C, "bo", linewidth=1.5)
        ax.plot(0, self.material.sigma_y, "bo", linewidth=1.5)
        ax.plot(self.material.sigma_y, 0, "bo", linewidth=1.5)
        ax.plot(self.material.sigma_u, 0, "bo", linewidth=1.5)
        ax.plot(0, P_y, "bo", linewidth=1.5)
        ax.plot(np.array([0, P_x]), np.array([0, P_y]), "g--", linewidth=1.5)
        ax.plot(np.array([0, P_x]), np.array([P_y, P_y]), "g--", linewidth=1.5)
        ax.text(self.sigma_m_eq, self.sigma_a_eq - 35, "P" + self.name[1])
        ax.text(10, self.material.sigma_Dm1C + 15, r"$\sigma_{D-1}^{C}$")
        ax.text(self.material.sigma_y, 20, r"$\sigma_{y}$")
        ax.text(10, self.material.sigma_y + 15, r"$\sigma_{y}$")
        ax.text(self.material.sigma_u, 20, r"$\sigma_{u}$")
        ax.text(10, P_y - 35, r"$\sigma_{D, lim}$")
        ax.set_xlabel(r"Mean Stress - $\sigma_{m, eq}$ [MPa]")
        ax.set_ylabel(r"Alternating Stress - $\sigma_{a, eq}$ [MPa]")
        plt.title("Haigh Diagram @ " + self.name)
        plt.grid()
        plt.show()
    
    # Calculate Section Static Safety Factor
    def calculateSectionStaticSafetyFactor(self, profile=None):
        zV = self.loc
        for i in range(len(profile.locs)):
            if zV >= profile.locs[i] and zV < profile.locs[i + 1]:
                self.staticSF = profile.shaft.material.sigma_y / profile.sigma_id[i]
    
    # Calculate Section Mean and Alternating Stresses
    def calculateSectionMeanAlternatingStress(self, profile=None):
        zV = self.loc
        for j in range(len(profile.locs)):
            if zV >= profile.locs[j] and zV < profile.locs[j + 1]:
                self.sigma_m_N = profile.sigma_N[j]                    
                self.sigma_a_N = 0                    
                self.sigma_m_Mb = 0                    
                self.sigma_a_Mb = profile.sigma_Mb[j]                    
                self.tau_m_Mt = profile.tau_Mt[j]
                self.tau_a_Mt = 0
    
    # Calculate Fatigue Safety Factor
    def calculateSectionFatigueSafetyFactor(self):
        # Calculate coordinates of specific points
        coeff1 = np.polyfit(np.array([0, self.material.sigma_y]), np.array([self.material.sigma_y, 0]), deg=1)
        coeff2 = np.polyfit(np.array([0, self.material.sigma_u]), np.array([self.material.sigma_Dm1C, 0]), deg=1)
        inter_x = (coeff2[1] - coeff1[1]) / (coeff1[0] - coeff2[0])
        coeff3 = np.polyfit(np.array([0, self.sigma_m_eq]), np.array([0, self.sigma_a_eq]), deg=1)
        P_x = (coeff2[1] - coeff3[1]) / (coeff3[0] - coeff2[0])
        if P_x > inter_x:
            P_x = (coeff1[1] - coeff3[1]) / (coeff3[0] - coeff1[0])
        P_y = np.polyval(coeff3, P_x)
        self.fatigueSF = P_y / self.sigma_a_eq