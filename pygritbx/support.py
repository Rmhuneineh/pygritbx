'''
This is the "Support" class. It inherits "Component" class.
It defines the properties for SKF bearing supports:
1) "name": a string of characters acting as a label
2) "type": "Pin" / "Roller"
3) "catalogueName": name of bearing as specified in catalogue
4) "catalogueType": "Standard" / "Explorer"
5) "d": bearing internal diamater expressed in [mm]
6) "D": bearing external diameter expressed in [mm]
7) "dm": bearing mean diameter expressed in [mm]
8) "B": bearing width expressed in [mm]
9) "C": dynamic basic load rating expressed in [N]
10) "C0": static basic load rating expressed in [N]
11) "Pu": fatigue load limit expressed in [N]
12) "a": bearing pressure point offset expressed in [mm]
13) "e": calculation factor 1
14) "X": bearing radial load factor
15) "Y": bearing axial load factor
16) "Y0": calculation factor 2
17) "Y1": calculation factor 3
18) "Y2": calculation factor 4
19) "A": bearing calculation factor 1
20) "kr": minimium load factor
21) "shoulder": "1" for +ve offser / "-1" for negative offset
22) "Frm": minimium load expressed in [N]
23) "P": equivalent dynamic bearing load expressed in [N]
24) "P0": equivalent static bearing load expressed in [N]
25) "F_a": axial reaction force expressed in [N]
26) "F_r" radial reaction force expressed in [N]
27) "arr": "Single" / "F2F" / "B2B" / "Tandem" / "Double row"
28) "n": rotational speed expressed in [rpm]
29) "nr": reference speed expressed in [rpm]
30) "s0": safety factor for static bearing load
31) "p": exponent for life equation
32) "a1": reliability factor
33) "eta_c": contamination factor
34) "a_skf": lubrication, contamination, and fatigue limit factor
35) "L_10m": bearing life expressed in [millions of cycles]
36) "L_10mh": bearing life expressed in [hours]
'''
from .component import Component
from .force import Force
import numpy as np
class Support(Component):

    # Reference relibility factor
    rel_ref = np.array([90, 95, 96, 97, 98, 99])
    a1_ref = np.array([1, 0.64, 0.55, 0.47, 0.37, 0.25])

    # Constructor
    def __init__(self, name="", type="", bearingType="", catalogueName="", catalogueType="", d=0, D=0, B=0,
                C=0, C0=0, Pu=0, nr=0, a=0, e=0, X=0, Y=0, Y0=0, Y1=0, Y2=0, A=0, kr=0, shoulder=0, arr="",
                axis=np.zeros(3), loc=0):
        super().__init__(name=name, material=None, axis=axis, loc=loc)
        self.type = type
        self.bearingType = bearingType
        self.catalogueName = catalogueName
        self.catalogueType = catalogueType
        self.d = d
        self.D = D
        self.dm = (self.d + self.D) / 2
        self.B = B
        self.C = C
        self.C0 = C0
        self.Pu = Pu
        self.a = a
        self.e = e
        self.X = X
        self.Y = Y
        self.Y0 = Y0
        self.Y1 = Y1
        self.Y2 = Y2
        self.A = A
        self.kr = kr
        self.shoulder = shoulder
        self.nr = nr
        self.n = 0
        self.arr = arr
        if self.bearingType == "ball":
            self.p = 3
        else:
            self.p = 10 / 3
        self.a_skf = 0

    # Life analysis
    def performLifeAnalysis(self, rel=100, condition="", a_skf=0):
        self.a_skf = a_skf
        print(f"Initiating Life Analysis on bearing {self.name}.")
        print(f"Checking minimum load condition.")
        if self.calculateMinimumLoad():
            print(f"Calculating static safety factor.")
            self.calculateEquivalentStaticLoad()
            print(f"Bearing {self.name}'s equivalent static load: {self.P0:.2f} [N].")
            print(f"Bearing {self.name}'s static safety factor: {self.s0:.2f} [-].")
            print(f"Calculating reliability factor.")
            self.calculateA1(rel=rel)
            print(f"Bearing {self.name}'s reliability factor: {self.a1:.2f} [-].")
            print(f"Calculating equivalent dynamic load.")
            self.calculateEquivalentDynamicLoad()
            print(f"Bearing {self.name}'s equivalent dynamic load: {self.P:.2f} [N].")
            print(f"Calculating contamination factor based on given condition: '{condition}'.")
            self.calculateEtaC(condition=condition)
            print(f"Bearing {self.name}'s contamination factor: {self.eta_c:.2f} [-].")
            print(f"Calculating bearing life.")
            self.calculateBearingLife()
            print(f"Bearing {self.name} life analysis results: {self.L_10m:.2f} [million cycles] | {self.L_10mh:.2f} [hours].")
    
    # Update reaction force
    def updateReaction(self):
        F_aV = self.F_tot.force * self.axis
        self.F_a = np.sqrt(np.sum(F_aV * F_aV))
        F_rV = self.F_tot.force - self.F_a
        self.F_r = np.sqrt(np.sum(F_rV * F_rV))
    
    # Minimum load calculation
    def calculateMinimumLoad(self):
        if self.bearingType == "Tapered":
            if self.catalogueType == "Standard":
                self.Frm = 0.02 * self.C
            elif self.catalogueType == "Explorer":
                self.Frm = 0.017 * self.C
        elif self.bearingType == "Cylindrical":
            self.Frm = self.kr * (6 + 4 * abs(self.n) / self.nr) * (self.dm / 100) ** 2 * 1e3
        else:
            raise ValueError(f"Bearing type '{self.bearingType}' not available.")
        if self.F_r != None:
            if self.F_r >= self.Frm:
                print(f"Bearing {self.name} satisfies minimium load condition.")
                return True
            else:
                print(f"Bearing {self.name} does not satisfy minimium load condition.")
                return False
    
    # Calculate equivalent dynamic load
    def calculateEquivalentDynamicLoad(self):
        ratioCond = self.F_a / self.F_r <= self.e
        if self.bearingType == "Tapered":
            if self.arr == "Single":
                if ratioCond:
                    self.P = self.F_r
                else:
                    self.P = 0.4 * self.F_r + self.Y * self.F_a
            elif self.arr == "F2F":
                if ratioCond:
                    self.P = self.F_r + self.Y1 * self.F_a
                else:
                    self.P = 0.67 * self.F_r + self.Y2 * self.F_a
            elif self.arr == "B2B":
                if ratioCond:
                    self.P = self.F_r + self.Y1 * self.F_a
                else:
                    self.P = 0.67 * self.F_r + self.Y2 * self.F_a
            elif self.arr == "Tandem":
                if ratioCond:
                    self.P = self.F_r
                else:
                    self.P = 0.4 * self.F_r + self.Y * self.F_a
            elif self.arr == "Double row":
                if ratioCond:
                    self.P = self.F_r + self.Y1 * self.F_a
                else:
                    self.P = 0.67 * self.F_r + self.Y2 * self.F_a
            else:
                raise ValueError("Bearing arrangement not available.")
        elif self.bearingType == "Cylindrical":
            if ratioCond:
                self.P = self.F_r
            else:
                self.P = 0.92 * self.F_r + self.Y * self.F_a
        else:
            raise ValueError("Bearing type not available.")
    
    # Calculate equivalent static load
    def calculateEquivalentStaticLoad(self):
        if self.bearingType == "Tapered":
            if self.arr == "Single":
                self.P0 = 0.5 * self.F_r + self.Y0 * self.F_a
                if self.P0 < self.F_r:
                    self.P0 = self.F_r
            elif self.arr == "F2F":
                self.P0 = self.F_r + self.Y0 * self.F_a
                if self.P0 < self.F_r:
                    self.P0 = self.F_r
            elif self.arr == "B2B":
                self.P0 = self.F_r + self.Y0 * self.F_a
                if self.P0 < self.F_r:
                    self.P0 = self.F_r
            elif self.arr == "Tandem":
                self.P0 = 0.5 * self.F_r + self.Y0 * self.F_a
            elif self.arr == "Double row":
                self.P0 = self.F_r + self.Y0 * self.F_a
                if self.P0 < self.F_r:
                    self.P0 = self.F_r
            else:
                raise ValueError("Bearing arrangement not avaialabe.")
        if self.bearingType == "Cylindrical":
            self.P0 = self.F_r
        self.s0 = self.C0 / self.P0
    
    # Calculate a1
    def calculateA1(self, rel=100):
        self.a1 = np.interp(rel, self.__class__.rel_ref, self.__class__.a1_ref)
    
    # Calculate contamination factor
    def calculateEtaC(self, condition=""):
        if condition == "Extreme cleanliness":
            self.eta_c = 1
        elif condition == "High cleanliness":
            if self.dm < 100:
                self.eta_c = 0.7
            else:
                self.eta_c = 0.85
        elif condition == "Normal cleanliness":
            if self.dm < 100:
                self.eta_c = 0.55
            else:
                self.eta_c = 0.7
        elif condition == "Slight contamination":
            if self.dm < 100:
                self.eta_c = 0.4
            else:
                self.eta_c = 0.5
        elif condition == "Typical contamination":
            if self.dm < 100:
                self.eta_c = 0.2
            else:
                self.eta_c = 0.3
        elif condition == "Severe contamination":
            self.eta_c = 0.05
        elif condition == "Very severe contamination":
            self.eta_c = 0
        else:
            raise ValueError("Invalid contamination condition.")

    # Calculate bearing life
    def calculateBearingLife(self):
        self.L_10m = self.a1 * self.a_skf * (self.C / self.P) ** self.p
        self.L_10mh = 1e6 / 60 / np.abs(self.n) * self.L_10m