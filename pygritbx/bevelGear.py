import numpy as np
from math import pi, ceil, cos, sin
from .gear import Gear

class BevelGear(Gear):

    # Constructor
    def __init__(self, name="", axis=np.zeros(3), loc=0, m_n=0, z=0, gamma=0, phi_n=0, Q_v=0, FW=0, material=None):
        super().__init__(name=name, axis=axis, loc=loc, m_n=m_n, z=z, psi=0, phi_n=phi_n, Q_v=Q_v, FW=FW, material=material)
        # Given parameters
        self.gamma = gamma * pi / 180
        # Calculated parameters
        self.z_p = ceil(self.z / cos(self.gamma))
        self.d_av = self.d - self.FW * sin(np.abs(self.gamma))