import numpy as np
from math import pi
class ShaftProfile:

    # Constructor
    def __init__(self, name="", radii=np.array([]), locs=np.array([])):
        self.name = name
        # Geometry
        self.radii = np.concatenate(([0], radii, [0]))
        self.locs = np.concatenate(([locs[0]], locs, [locs[-1]]))
        # Shaft
        self.shaft = None
        # Cross-sectional properties
        self.Area = np.array([])
        self.Wb = np.array([])
        self.Wt = np.array([])
        # Stresses
        self.sigma_N = np.array([])
        self.sigma_Mb = np.array([])
        self.tau_Mt = np.array([])
        self.sigma_tot = np.array([])
        self.sigma_id = np.array([])
    
    # Overload Equal
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return all(self.radii == other.radii) and all(self.locs == other.locs)
        return False
    
    # Add Fillet
    def addFillet(self, radius=0, quadrant=[], zOff=0, dOff=0):
        for q in range(len(quadrant)):
            if quadrant[q] == 1:
                theta = np.arange(pi, 3*pi/2, 0.1)
            elif quadrant[q] == 2:
                theta = np.arange(3*pi/2, 2*pi, 0.1)
            else:
                raise ValueError("Wrong Input")
        z = radius * np.cos(theta) + zOff
        r = radius * np.sin(theta) + dOff
        before = self.locs < np.min(z)
        after = self.locs > np.max(z)
        self.radii = np.concatenate((self.radii[np.where(before)], r, self.radii[np.where(after)]))
        self.locs = np.concatenate((self.locs[np.where(before)], z, self.locs[np.where(after)]))
    
    # Refine Profile
    def refineProfile(self, name="", delta=0.1):
        if name == "":
            name = self.name
        pLen = int((self.locs[-1] - self.locs[0]) / delta + 1)
        refinedProfile = ShaftProfile(name, np.zeros(pLen), np.arange(self.locs[0], self.locs[-1] + delta / 2, delta))
        for z in range(2, len(self.locs) - 3):
            condition = np.where(np.logical_and(refinedProfile.locs >= self.locs[z], refinedProfile.locs <= self.locs[z + 1]))
            if self.locs[z] != self.locs[z + 1]:
                refinedProfile.radii[condition] = np.interp(refinedProfile.locs[condition], np.array([self.locs[z], self.locs[z+1]]), np.array([self.radii[z], self.radii[z+1]]))
        return refinedProfile
    
    # Calculate Cross-Sectional Properties
    def calculateSectionProperties(self):
        self.Area = pi * self.radii ** 2
        self.Wb = pi / 4 * self.radii ** 3
        self.Wt = pi / 2 * self.radii ** 3
    
    # Plot Profile
    def plotProfile(self, ax=None):
        ax1 = ax.twinx()
        ax1.plot(self.locs, self.radii, 'r', linewidth=1.5)
        ax1.plot(self.locs, -self.radii, 'r', linewidth=1.5)
        window = (self.locs[-1] - self.locs[0] + 20) / 2
        ax1.set_xlim(-0.1 * window, 2 * window)
        ax1.set_ylim(-window, window)
        ax1.set_ylabel("Profile [mm]")
    
    # Calculate Profile Stresses
    def calculateProfileStresses(self):
        if self.sigma_N.size != 0 and self.sigma_Mb.size != 0 and self.tau_Mt.size != 0:
            return
        sLen = len(self.locs)
        self.sigma_N = np.zeros(sLen)
        self.sigma_Mb = np.zeros(sLen)
        self.tau_Mt = np.zeros(sLen)
        self.sigma_N[np.where(self.Area != 0)] = self.shaft.N[np.where(self.Area != 0)] / self.Area[np.where(self.Area != 0)]
        self.sigma_Mb[np.where(self.Wb != 0)] = 1e3 * self.shaft.Mf[np.where(self.Wb != 0)] / self.Wb[np.where(self.Wb != 0)]
        self.tau_Mt[np.where(self.Wt != 0)] = 1e3 * self.shaft.Mt[np.where(self.Wt != 0)] / self.Wt[np.where(self.Wt != 0)]
    
    # Calculate Profile Equivalent and Ideal Stresses
    def calculateProfileEquivalentAndIdealStress(self):
        if self.sigma_tot.size != 0 and self.sigma_id.size != 0:
            return
        self.sigma_tot = self.sigma_N + self.sigma_Mb
        self.sigma_id = np.sqrt(self.sigma_tot ** 2 + 3 * self.tau_Mt ** 2)
    
    # Plot stresses
    def plotStresses(self):
        # Normal stress
        self.shaft.plotLoad(load=self.sigma_N, ylabel=r"$\sigma^{N}$ [MPa]", title=r"Normal Stress - $\sigma^{N}(z)$ [MPa]", profile=self)
        # Bending stress
        self.shaft.plotLoad(load=self.sigma_Mb, ylabel=r"$\sigma^{M_{B}}$ [MPa]", title=r"Bending Stress - $\sigma^{M_{B}}(z)$ [MPa]", profile=self)
        # Torsional stress
        self.shaft.plotLoad(load=self.tau_Mt, ylabel=r"$\tau^{M_{t}}$ [MPa]", title=r"Torsional Stress - $\tau^{M_{t}}(z)$ [MPa]", profile=self)
        # Total stress
        self.shaft.plotLoad(load=self.sigma_tot, ylabel=r"$\sigma^{tot}$ [MPa]", title=r"Resulting Normal Stress - $\sigma^{tot}(z)$ [MPa]", profile=self)
        # Equivalent stress
        self.shaft.plotLoad(load=self.sigma_id, ylabel=r"$\sigma_{id}$ [MPa]", title=r"Equivalent Stress - $\sigma_{id}(z)$ [MPa]", profile=self)