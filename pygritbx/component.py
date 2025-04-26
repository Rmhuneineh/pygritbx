import numpy as np
'''
This is the "Component" class.
It's a parent class that defines general properties that are common among different components.
Other classes will inherit these properties instead of having to define every single time.
Properties are:
1) "name": a string of characters acting as a label
2) "material": a "Material" object of the material class defining the material properties of the component
3) "axis": a 3-element vector representing the axis along which the component is rotating with respect to a defined reference frame
4) "loc": a 3-element vector representing the location of the component with respect to a defined reference frame
5) "F_tot": the total force acting on the component expressed in [N]
6) "T_tot": the total torque acting on the component expressed in [Nm]
'''
class Component:

    # Constructor
    def __init__(self, name, material, axis, loc, EFs=np.array([]), ETs=np.array([]), omega=np.zeros(3)):
        self.name = name
        self.material = material
        self.axis = axis
        self.loc = loc
        self.EFs = EFs
        self.ETs = ETs
        self.omega = omega
    
    # Check force equilibrium
    def checkForceEquilibrium(self):
        eq = np.zeros(3)
        for EF in self.EFs:
            eq = eq + EF.force
        if all(eq <= 1e-3 * np.ones(3)):
            print(f"{self.name} maintains a force equilibrium.")
        else:
            print(f"{self.name} does not maintain a force equilibrium.")
    
    # Check torque equilibrium
    def checkTorqueEquilibrium(self):
        pass
    
    # Update external forces
    def updateEFs(self, EFs):
        for ef in EFs:
            if ef not in self.EFs:
                self.EFs = np.append(self.EFs, ef)
    
    # Update external torques
    def updateETs(self, ETs):
        for et in ETs:
            if et not in self.ETs:
                self.ETs = np.append(self.ETs, et)