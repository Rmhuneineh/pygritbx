'''
This is the "Motor" class. It inherits from "Component" class.
The motor class defines the properties of an electric motor based on:

I) Given properties:
--> 1) "name": a string of characters acting as a label
--> 2) "power": scalar value representing the power of the motor expressed in [W]
--> 3) "n": scalar value representing the speed expressed in [rpm]
--> 4) "torque": a 3-element vector representing the torque exerted by the motor expressed in [Nm]
 For properties (2), (3), and (4), the user can indicate either 2 of them (the 3rd will be calculated) or all of them where a check will be performed.

II) Calculated properties
--> 1) "omega": a 3-element vector representing the rotational velocity of the motor expressed in [rad/s]

The functions implemented help relocate the load(s) exerted by the motor once placed on a shaft.
'''
import numpy as np
from .component import Component
from .force import Force
from .torque import Torque
from math import pi

class Motor(Component):
    
    # Constructor
    def __init__(self, name="", loc=0, power=0, n=0, torque=Torque(np.array([]), np.array([])), axis=np.array([0, 0, 0])):
        # Given properties
        super().__init__(name=name, material=None, axis=axis, loc=loc)
        # Check for valid input
        if power == 0 and n == 0 and torque.size == 0:
            raise ValueError("A minimum of 2 out of 3 inputs are necessary: power, n, and torque.")
        elif power != 0 and n != 0 and torque.size != 0:
            omega = n * pi / 30 * self.axis
            omega_mag = np.sqrt(np.sum(omega * omega))
            if power != torque.mag() * omega_mag:
                raise ValueError(f"Provided given is incoherent: {power} != {torque.mag() * omega_mag}.")
        # Check which inputs are given
        if power == 0:
            omega = n * pi / 30 * self.axis
            omega_mag = np.sqrt(np.sum(omega * omega))
            power = torque.mag() * omega_mag
        elif n == 0:
            omega = power / torque.mag() * self.axis
            omega_mag = np.sqrt(np.sum(omega * omega))
            n = omega_mag * 30 / pi
        elif torque.torque.size == 0:
            omega = n * pi / 30
            torque_mag = power / omega
            torque = torque_mag * self.axis

        self.power = power
        self.n = n
        self.omega = omega
        if self.abs_loc.size != 0:
            location = self.abs_loc
        else:
            location = self.rel_loc
        self.updateETs([Torque(torque=torque, loc=location)])
        self.updateEFs([Force(force=np.array([0, 0, 0]), loc=location)])
    
    # Update Force Location
    def updateForceLoc(self):
        for EF in self.EFs:
            EF.loc = self.abs_loc
    
    # Update Torque Location
    def updateTorqueLoc(self):
        for ET in self.ETs:
            ET.loc = self.abs_loc