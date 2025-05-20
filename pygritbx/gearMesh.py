'''
This is the "Gear Mesh" class.
It defines a gear mesh by relying on 2 gear components and the type of mesh.

The properties used are:
I) Given properties
--> 1) "drivingGear": a gear component representing the input gear of the gear mesh
--> 2) "drivenGear": a gear component represnting the output gear of the gear mesh
--> 3) "axis": a 3-element vector representing the direction along which the centers of the driving gear and the driven gear centers are connected
--> 4) "type": a string of characters indicating the gear mesh type: Internal / External

II) Calculated properties
--> 1) "ratio": the gear ratio calculated based on the number of teeth
--> 2) "m_G": the reciprocal of the gear ratio
--> 3) "loc": a 3-element vector representing the location of the gear mesh with respect to the origin [0, 0, 0]
--> 4) "F": a 3-element vector representing the resultant force acting on the driven gear due to the gear mesh object defined
--> 5) "F_t": a 3-element vector representing the tangential force component of the resultant force
--> 6) "F_r": a 3-element vector representing the radial force component of the resultant force
--> 7) "F_a": a 3-element vector representing the axial force component of the resultant force

'''
import numpy as np
from .torque import Torque
from .force import Force
class GearMesh:

    # Constructor
    def __init__(self, name="", drivingGear=None, drivenGear=None, axis=np.zeros(3), type=""):
        if drivingGear == None or drivenGear == None:
            raise ValueError("Driving or driven gear missing!")
        if drivingGear.m_n != drivenGear.m_n:
            raise Exception("Incompatible Gear Mesh!")
        # Update location of driven gear
        if drivenGear.abs_loc.size == 0:
            drivenGear.abs_loc = drivingGear.abs_loc + axis * (drivingGear.d + drivenGear.d) / 2
        # Given properties
        self.name = name
        self.drivingGear = drivingGear
        self.drivenGear = drivenGear
        self.axis = axis
        self.type = type
        # Calculated properties
        self.ratio = self.drivingGear.z / self.drivenGear.z
        self.m_G = 1 / self.ratio
        sgn = 1 # aassuming self.type = "External"
        if self.type == "External":
            sgn = -1
        self.drivenGear.omega = sgn * self.ratio * self.drivingGear.omega
        #self.drivenGear.T_tot = Torque(-sgn * self.drivingGear.T_tot.torque / self.ratio, self.drivenGear.loc)
        self.loc = self.drivingGear.d / 2 * self.axis + self.drivingGear.abs_loc
        self.F = Force(np.zeros(3), self.loc) # Resultant Force
        self.F_t = Force(np.zeros(3), self.loc) # Tangential Force
        self.F_r = Force(np.zeros(3), self.loc) # Radial Force
        self.F_a = Force(np.zeros(3), self.loc) # Axial Force
        # Update gear meshes
        self.drivingGear.meshes = np.append(self.drivingGear.meshes, self)
        self.drivenGear.meshes = np.append(self.drivenGear.meshes, self)