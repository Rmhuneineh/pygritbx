import numpy as np
from component import Component
from force import Force
from torque import Torque
from math import pi

class InputMotor(Component):

    def __init__(self, name, power, n, axis, loc):
        super().__init__(name=name, material=None, axis=axis, loc=loc, F_tot=None, T_tot=None, omega=None)
        self.power = power
        self.n = n
        self.omega = self.n * pi / 30 * self.axis
        self.T_tot = Torque(np.array([0 if o == 0 else self.power/o for o in self.omega]), self.loc)
        self.F_tot = Force(np.array([0, 0, 0]), self.loc)
        
