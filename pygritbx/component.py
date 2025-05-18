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
    def __init__(self, name="", material=None, axis=np.zeros(3), loc=0, EFs=np.array([]), ETs=np.array([]), omega=np.zeros(3)):
        self.name = name
        self.material = material
        self.axis = axis
        if not isinstance(loc, list):
            self.rel_loc = loc * self.axis
            self.abs_loc = np.array([])
        else:
            self.rel_loc = np.array([])
            self.abs_loc = np.array(loc)
        self.EFs = EFs
        self.ETs = ETs
        self.omega = omega
        self.onShaft = None
    
    # Check force equilibrium
    def checkForceEquilibrium(self):
        print(f"Checking force equilibrium on {self.name}.")
        valid = False
        if self.EFs.size != 0:
            valid = True
        else:
            return valid
        eq = np.zeros(3)
        eqState = False
        for EF in self.EFs:
            eq += EF.force
        if all(np.abs(eq) <= 1e-3 * np.ones(3)):
            print(f"{self.name} maintains a force equilibrium.")
            eqState = True
        else:
            print(f"{self.name} does not maintain a force equilibrium.")
        return eqState

    # Check torque equilibrium
    def checkTorqueEquilibrium(self):
        pass

    # Solve function
    def solve(self):
        pass
    
    # Update external forces
    def updateEFs(self, EFs=[]):
        for ef in EFs:
            if ef not in self.EFs:
                self.EFs = np.append(self.EFs, ef)
    
    # Update external torques
    def updateETs(self, ETs=[]):
        for et in ETs:
            if et not in self.ETs:
                self.ETs = np.append(self.ETs, et)
    
    # Update location
    def updateLoc(self):
        if self.abs_loc.size == 0:
            self.abs_loc = self.rel_loc + self.onShaft.abs_loc
        elif self.rel_loc.size == 0:
            self.rel_loc = self.abs_loc - self.onShaft.abs_loc
