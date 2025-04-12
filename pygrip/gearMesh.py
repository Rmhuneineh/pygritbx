'''
This is the gearMesh class.
It defines a gear mesh by relying on 2 gear components and the type of mesh.
The properties used are:
I) Given properties
--> 1) "drivingGear": a gear component to which the input is connected to
--> 2) "drivenGear": a gear component to which the output in connectec to
--> 3) "axis": a 3-element vector representing the axis of rotation
--> 4) "type": a string of characters indicating the gear mesh type: Internal / External
II) Calculated properties
--> 1) "ratio": the gear ratio calculated based on the number of teeth
--> 2) "m_G": the reciprocal of the gear ratio
After creating an instance of the class, it's possible to get the rotational velocity as well as the torque corresponding to the driven gear via:
1) GetOmegaMesh()
2) GetMeshTorque()
methods, respectively.
'''
from .torque import Torque
class GearMesh:

    # Constructor
    def __init__(self, drivingGear, drivenGear, axis, type):
        if drivingGear.m_n != drivenGear.m_n:
            raise Exception("Incompatible Gear Mesh!")
        # Given properties
        self.drivingGear = drivingGear
        self.drivenGear = drivenGear
        self.axis = axis
        self.type = type
        # Calculated properties
        self.ratio = self.drivingGear.z / self.drivenGear.z
        self.m_G = 1 / self.ratio
        sgn = -1 # aassuming self.type = "External"
        if self.type == "Internal":
            sgn = 1
        self.drivenGear.omega = sgn * self.ratio * self.drivingGear.omega
        self.drivenGear.T_tot = Torque(-sgn * self.drivingGear.T_tot.torque / self.ratio, self.drivenGear.loc)
    
    # Get Driven Gear Omega
    def GetOmegaMesh(self):
        return self.drivenGear.omega
    
    # Get Driven Gear Torque
    def GetMeshTorque(self):
        return self.drivenGear.T_tot.torque