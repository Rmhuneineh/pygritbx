import numpy as np
import matplotlib.pyplot as plt
from .component import Component
from .torque import Torque
from .force import Force
from .gear import Gear
from .motor import Motor
from math import pi
class Shaft(Component):

    # Constructor
    def __init__(self, name, inputs, outputs, axis, material, sups, loc=0):
        # Given parameters
        super().__init__(name=name, material=material, axis=axis, loc=loc, omega=inputs[0].omega)
        # Update shaft's absolute location
        if self.abs_loc.size == 0:
            self.abs_loc = -inputs[0].rel_loc + inputs[0].abs_loc
        # Inputs
        for input in inputs:
            input.onShaft = self
            input.updateLoc()
            if isinstance(input, Motor):
                input.updateForceLoc()
                input.updateTorqueLoc()
            if input.ETs.size != 0:
                self.updateETs(input.ETs)
        self.inputs = inputs
        # Outputs
        for out in outputs:
            out.onShaft = self
            out.updateLoc()
            if out.ETs.size != 0:
                self.updateETs(out.ETs)
            out.omega = self.omega
        self.outputs = outputs
        
        # Supports
        for sup in sups:
            sup.onShaft = self
            sup.updateLoc()
            sup.omega = self.omega
            sup.n = sup.omega * 30 / pi
        self.supports = sups
        # Sections
        self.sections = np.array([])
    
    # Solve function
    def solve(self):
        unknown_Ts = 0
        unknkown_Fs = 0
        unknown_comp_T = None
        comps = np.append(self.inputs, self.outputs)
        if not self.checkTorqueEquilibrium():
            print(f"Checking solvability for {self.name}.")
            for comp in comps:
                if comp.ETs.size == 0:
                    unknown_Ts += 1
                    unknown_comp_T = comp
                else:
                    self.updateETs(comp.ETs)
                if comp.EFs.size == 0 and not isinstance(comp, Motor):
                    unknkown_Fs += 1
                else:
                    self.updateEFs(comp.EFs)
            if unknown_Ts > 1:
                print(f"{self.name}'s torque equilibrium cannot be solved.")
            elif unknown_Ts == 1:
                print(f"Solving torque equilibrium for {self.name}.")
                self.calculateTorque(unknown_comp_T)
                print(f"Torque equilibrium for {self.name} is solved.")
                if isinstance(unknown_comp_T, Gear):
                    choice = input(f"Detected torque on {unknown_comp_T.name}. Would you like to solve its torque equilibrium [y/n]: ")
                    if choice == 'y' or choice == 'Y':
                        unknown_comp_T.solve()
                        unknkown_Fs -= 1
                    else:
                        print(f"{unknown_comp_T.name}'s torque equilibrium won't be solved now.")
        else:
            print(f"No torque equilibrium to be solved for {self.name}.")
        if unknkown_Fs == 0:
            if self.checkForceEquilibrium() :
                print(f"No force equilibrium to be solved for {self.name}.")
            else:
                reaction_choice = input(f"Forces from external components on {self.name} are resolved. Do you want to calculate the reaction forces [y/n]: ")
                if reaction_choice == 'y' or reaction_choice == 'Y':
                    self.calculateReactionForces()
                    self.checkForceEquilibrium()
                else:
                    print(f"Reaction forces on {self.name} won't be calculated.")
        else:
            print(f"Force equilibrium for {self.name} cannot be solved.")
            
        
    # Check torque equilibrium
    def checkTorqueEquilibrium(self):
        print(f"Checking torque equilibrium for {self.name}.")
        valid = True
        if self.ETs.size == 0:
            valid = False
            return valid
        eq = np.zeros(3)
        eqState = False
        for ET in self.ETs:
            eq = eq + ET.torque
        if all(np.abs(eq) <= 1e-3 * np.ones(3)):
            print(f"{self.name} mainatains a torque equilibrium.")
            eqState = True
        else:
            print(f"{self.name} does not mainatain a torque equilibrium.")
        return eqState
    
    # Calculate torque
    def calculateTorque(self, comp):
        ET = Torque(-np.sum(self.ETs), comp.abs_loc)        
        comp.updateETs([ET])
        self.updateETs([ET])
    
    # Set shaft profile
    def setProfile(self, profile):
        self.profile = profile
    
    # Add sections
    def addSections(self, sections):
        for section in sections:
            section.material = self.material
        self.sections = np.append(self.sections, sections)
    
    # Insert fatigue limit corrector factors
    def insertFLCF(self, sections):
        for sec1, sec2 in zip(self.sections, sections):
            sec1.AddFLCF()
            sec2.AddFLCF()
    
    # Calculate reaction forces
    def calculateReactionForces(self):
        # Gears axial load
        K_a = np.zeros(3)
        for EF in self.EFs:
            K_a = K_a + (EF.force * self.axis)
        # Find the bearing around which to apply moment
        index = 0
        for i in range(len(self.supports)):
            if self.supports[i].type == "Pin":
                index = i
                self.supports[index].F_tot = Force(np.zeros(3), self.supports[index].abs_loc)
        # Calculate reaction on other bearing
        for i in range(len(self.supports)):
            if i != index:
                self.supports[i].F_tot = Force(np.zeros(3), self.supports[i].abs_loc)
                supDist = self.supports[i].abs_loc - self.supports[index].abs_loc
                supDist_rec = np.array([1/d if d != 0 else 0 for d in supDist])
                for EF in self.EFs:
                    momEF = np.cross(EF.force, (EF.loc - self.supports[index].abs_loc)) * 1e-3
                    self.supports[i].F_tot.force += np.cross(momEF, supDist_rec) * 1e3
                self.updateEFs([self.supports[i].F_tot])
        # Calculate reaction around bearing with sum of external forces
        for EF in self.EFs:
            self.supports[index].F_tot.force -= EF.force
        self.updateEFs([self.supports[index].F_tot])
        # Update axial load based on configuration
        if self.supports[0].bearingType == "Tapered" and self.supports[1].bearingType == "Tapered":
            if self.supports[0].shoulder == -1:
                indA = 1
                indB = 0
                sgn = 1
            else:
                indA = 0
                indB = 1
                sgn = -1
            A_FrV = self.supports[indA].F_tot.force - self.supports[indA].F_tot.force * self.axis
            A_Fr = np.sqrt(np.sum(A_FrV * A_FrV))
            A_Y = self.supports[indA].Y
            B_FrV = self.supports[indB].F_tot.force - self.supports[indB].F_tot.force * self.axis
            B_Fr = np.sqrt(np.sum(B_FrV * B_FrV))
            B_Y = self.supports[indB].Y
            print(f"Axial reaction forces on {self.name}: ", end="")
            # Case 1
            if np.sum(K_a) > 0:
                # Factor of comparison
                fac = B_Fr / B_Y - A_Fr / A_Y
                # Case 1a
                if fac <= 0 and np.abs(np.sum(K_a)) >= 0:
                    print("Case 1a")
                    A_Fa = sgn * 0.5 * A_Fr / A_Y * self.axis
                    B_Fa = -(A_Fa + K_a)
                # Case 1b
                elif fac > 0 and np.abs(np.sum(K_a)) >= 0.5 * fac:
                    print("Case 1b")
                    A_Fa = sgn * 0.5 * A_Fr / A_Y * self.axis
                    B_Fa = -(A_Fa + K_a)
                # Case 1c
                elif fac > 0 and np.abs(np.sum(K_a)) < 0.5 * fac:
                    print("Case 1c")
                    B_Fa = -sgn * 0.5 * B_Fr / B_Y * self.axis
                    A_Fa = -(B_Fa + K_a)
            # Case 2
            else:
                # Factor of comparison
                fac = A_Fr / A_Y - B_Fr / B_Y
                # Case 2a
                if fac <= 0 and np.abs(np.sum(K_a)) >= 0:
                    print("Case 2a")
                    B_Fa = -sgn * 0.5 * B_Fr / B_Y * self.axis
                    A_Fa = -(B_Fa + K_a)
                # Case 2b
                elif fac > 0 and np.abs(np.sum(K_a)) >= 0.5 * fac:
                    print("Case 2b")
                    B_Fa = -sgn * 0.5 * B_Fr / B_Y * self.axis
                    A_Fa = -(B_Fa + K_a)
                # Case 2c
                elif fac > 0 and np.abs(np.sum(K_a)) < 0.5 * fac:
                    print("Case 2c")
                    A_Fa = sgn * 0.5 * A_Fr / A_Y * self.axis
                    B_Fa = -sgn*(A_Fa + K_a)
            self.supports[indA].F_tot.force[2] = np.sum(A_Fa)
            self.supports[indB].F_tot.force[2] = np.sum(B_Fa)
            if self.supports[0].arr == "B2B":
                self.EFs[-2].force[2] = np.sum(B_Fa)
                self.EFs[-1].force[2] = np.sum(A_Fa)
            elif self.supports[0].arr == "F2F":
                self.EFs[-2].force[2] = np.sum(A_Fa)
                self.EFs[-1].force[2] = np.sum(B_Fa)
        # Update support reaction to separate total radial force and axial force
        for support in self.supports:
            support.updateReaction()
    
    # Calculate internal loads
    def calculateInternalLoads(self, RF):
        l = len(self.profile.locs)
        self.N = np.zeros(l)
        self.Mx = np.zeros(l)
        self.My = np.zeros(l)
        self.Mt = np.zeros(l)
        for EF in self.EFs:
            for i, z in enumerate(self.profile.locs):
                if np.dot(EF.loc - self.abs_loc, np.abs(self.axis)) <= z:
                    self.N[i] = self.N[i] - np.sum(EF.force * np.abs(self.axis))
                    mxz = np.sum(np.cross(EF.force * RF[2], EF.loc * RF[1]))
                    mxy = np.sum(np.cross(EF.force * RF[1], (z - EF.loc ) * RF[2]))
                    self.Mx[i] = self.Mx[i] + (mxz - mxy) * 1e-3
                    myz = np.sum(np.cross(EF.force * RF[2], EF.loc * RF[0]))
                    myx = np.sum(np.cross(EF.force * RF[0], (EF.loc - z) * RF[2]))
                    self.My[i] = self.My[i] + (myz + myx) * 1e-3
        for ET in self.ETs:
                for i, z in enumerate(self.profile.locs):
                    if np.dot(ET.loc - self.abs_loc, np.abs(self.axis)) <= z:
                        self.Mt[i] = self.Mt[i] + np.sum(ET.torque)
        self.N[np.where(np.abs(self.N) < 1e-3)] = 0
        self.Mx[np.where(np.abs(self.Mx) < 1e-3)] = 0
        self.My[np.where(np.abs(self.My) < 1e-3)] = 0
        self.Mt[np.where(np.abs(self.Mt) < 1e-3)] = 0
        self.Mf = np.sqrt(self.Mx ** 2 + self.My ** 2)
    
    # Calculate stresses
    def calculateStresses(self):
        sLen = len(self.profile.locs)
        self.sigma_N = np.zeros(sLen)
        self.sigma_Mb = np.zeros(sLen)
        self.tau_Mt = np.zeros(sLen)
        self.sigma_N[np.where(self.profile.Area != 0)] = self.N[np.where(self.profile.Area != 0)] / self.profile.Area[np.where(self.profile.Area != 0)]
        self.sigma_Mb[np.where(self.profile.Wb != 0)] = 1e3 * self.Mf[np.where(self.profile.Wb != 0)] / self.profile.Wb[np.where(self.profile.Wb != 0)]
        self.tau_Mt[np.where(self.profile.Wt != 0)] = 1e3 * self.Mt[np.where(self.profile.Wt != 0)] / self.profile.Wt[np.where(self.profile.Wt != 0)]
    
    # Calculate equivalent and ideal stresses
    def calculateEquivalentAndIdealStress(self):
        self.sigma_tot = self.sigma_N + self.sigma_Mb
        self.sigma_id = np.sqrt(self.sigma_tot ** 2 + 3 * self.tau_Mt ** 2)
    
    # Plot internal loads
    def plotInternalLoads(self):
        # Normal load
        self.plotLoad(self.N, "N [N]", "Normal Load - N(z)")
        # Bending moment around x-axis
        self.plotLoad(self.Mx, r"$M_{x}$ [Nm]", r"Bending Moment $M_{x}(z)$")
        # Bending moment around y-axis
        self.plotLoad(self.My, r"$M_{y}$ [Nm]", r"Bending Moment $M_{y}(z)$")
        # Resulting bending moment
        self.plotLoad(self.Mf, r"$M_{B}$ [Nm]", r"Bending Moment $M_{B}(z)$")
        # Torsional moment
        self.plotLoad(self.Mt, r"$M_{t}$ [Nm]", r"Torsional Moment $M_{t}(z)$")
    
    # Plot load with shaft profile
    def plotLoad(self, load, ylabel, title):
        fig, ax = plt.subplots()
        ax.plot(self.profile.locs, load, 'b', linewidth = 1.5)
        ax.set_xlabel("z [mm]")
        ax.set_ylabel(ylabel)
        plt.title(title)
        if np.max(np.abs(load)) != 0:
            ax.set_ylim(-1.1 * np.max(np.abs(load)), 1.1 * np.max(np.abs(load)))
        else:
            ax.set_ylim(-1, 1)
        plt.grid()
        self.profile.plotProfile(ax)
        for section in self.sections:
            xs = np.ones(2) * np.sum(section.loc)
            ys = np.array([-1, 1])
            if np.max(np.abs(load)) != 0:
                ys = ys * 1.1 * np.max(np.abs(load))
            ax.plot(xs, ys, 'g--', linewidth=1.5)
            ax.text(xs[0] - 10, 0.9 * ys[0], section.name)
            ax.text(xs[0] - 10, 0.9 * ys[1], section.name)
    
    # Plot stresses
    def plotStresses(self):
        # Normal stress
        self.plotLoad(self.sigma_N, r"$\sigma^{N}$ [MPa]", r"Normal Stress - $\sigma^{N}(z)$ [MPa]")
        # Bending stress
        self.plotLoad(self.sigma_Mb, r"$\sigma^{M_{B}}$ [MPa]", r"Bending Stress - $\sigma^{M_{B}}(z)$ [MPa]")
        # Torsional stress
        self.plotLoad(self.tau_Mt, r"$\tau^{M_{t}}$ [MPa]", r"Torsional Stress - $\tau^{M_{t}}(z)$ [MPa]")
        # Total stress
        self.plotLoad(self.sigma_tot, r"$\sigma^{tot}$ [MPa]", r"Resulting Normal Stress - $\sigma^{tot}(z)$ [MPa]")
        # Equivalent stress
        self.plotLoad(self.sigma_id, r"$\sigma_{id}$ [MPa]", r"Equivalent Stress - $\sigma_{id}(z)$ [MPa]")
    
    # Calculate sections static safety factor
    def calculateStaticSafetyFactor(self, sections):
        for i in range(len(self.sections)):
            zV = np.sum(self.sections[i].loc)
            for j in range(len(self.profile.locs)):
                if zV >= self.profile.locs[j] and zV < self.profile.locs[j + 1]:
                    self.sections[i].staticSF = self.material.sigma_y / self.sigma_id[j]
            sections[i].staticSF = self.sections[i].staticSF

    # Calculate mean and alternating stresses
    def calculateMeanAlternatingStress(self, sections):
        for i in range(len(self.sections)):
            zV = np.sum(self.sections[i].loc)
            for j in range(len(self.profile.locs)):
                if zV >= self.profile.locs[j] and zV < self.profile.locs[j + 1]:
                    self.sections[i].sigma_m_N = self.sigma_N[j]
                    sections[i].sigma_m_N = self.sigma_N[j]
                    self.sections[i].sigma_a_N = 0
                    sections[i].sigma_a_N = 0
                    self.sections[i].sigma_m_Mb = 0
                    sections[i].sigma_m_Mb = 0
                    self.sections[i].sigma_a_Mb = self.sigma_Mb[j]
                    sections[i].sigma_a_Mb = self.sigma_Mb[j]
                    self.sections[i].tau_m_Mt = self.tau_Mt[j]
                    sections[i].tau_m_Mt = self.tau_Mt[j]
                    self.sections[i].tau_a_Mt = 0
                    sections[i].tau_a_Mt = 0
            sections[i].staticSF = self.sections[i].staticSF
    
    # Calculate equivalent mean and alternating stress
    def calculateEquivalentStresses(self, sections):
        for section in sections:
            section.CalculateEquivalentStress()