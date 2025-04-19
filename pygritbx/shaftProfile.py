import numpy as np
from math import pi
class ShaftProfile:

    # Constructor
    def __init__(self, radii, locs):
        self.radii = np.concatenate(([0], radii, [0]))
        self.locs = np.concatenate(([locs[0]], locs, [locs[-1]]))
    
    # Add Fillet
    def AddFillet(self, radius, quadrant, zOff, dOff):
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
    def refineProfile(self, delta):
        pLen = int((self.locs[-1] - self.locs[0]) / delta + 1)
        refinedProfile = ShaftProfile(np.zeros(pLen), np.arange(self.locs[0], self.locs[-1] + delta / 2, delta))
        for z in range(2, len(self.locs) - 3):
            condition = np.where(np.logical_and(refinedProfile.locs >= self.locs[z], refinedProfile.locs <= self.locs[z + 1]))
            if self.locs[z] != self.locs[z + 1]:
                refinedProfile.radii[condition] = np.interp(refinedProfile.locs[condition], np.array([self.locs[z], self.locs[z+1]]), np.array([self.radii[z], self.radii[z+1]]))
        return refinedProfile
    
    # Calculate Cross-Sectional Properties
    def CalculateSectionProperties(self):
        self.Area = pi * self.radii ** 2
        self.Wb = pi / 4 * self.radii ** 3
        self.Wt = pi / 2 * self.radii ** 3
    
    # Plot Profile
    def plotProfile(self, ax):
        ax1 = ax.twinx()
        ax1.plot(self.locs, self.radii, 'r', linewidth=1.5)
        ax1.plot(self.locs, -self.radii, 'r', linewidth=1.5)
        window = (self.locs[-1] - self.locs[0] + 20) / 2
        ax1.set_xlim(-0.1 * window, 2 * window)
        ax1.set_ylim(-window, window)
        ax1.set_ylabel("Profile [mm]")