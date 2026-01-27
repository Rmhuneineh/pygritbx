# Python-based Gearbox Reliability and Integrity Toolbox (PyGRITbx)
PyGritbx is a python-based tool born thanks to a project in the course **"Fundamentals of Machine Design"** at **Politecnico di Torino**.
Students following 3rd year of the Mechanical Engineering bachelor's degree are required to submit a report for this project.
The tool helps students define the components that constitute the gearbox in a given project, in order to then verify its design.

## What does it do?
Given a gearbox configuration and an operating point, the tool is able to:
1) Calculate all the forces exchanged between gears and reaction forces produced by the bearings.
2) Calculate internal loads and stresses on shafts.
3) Perform static and fatigue verification on shafts by calculating the corresponding safety factors.
4) Gear tooth verification in terms of bending and pitting (wear).
5) SKF bearing life analysis.

## How does it do it?
Thanks to **Python Object-Oriented Programming (POOP)**, the tool is able to define the components of the gearbox as objects as well as account for their specific geometries and interactions to perform the necessary calculations.

### Components
The tool supports the following components:
1) Motor
2) Gear (can be used as spur/helical gear by setting the proper helix angle) *and recently, even Bevel gears*
3) Shaft
4) SKF Bearings

This structure allows the user to define the components based on their given characteristics.

### Interactions and Geometries
The tool supports other classes that the user needs in order to define the interactions between different components as well as their geometry:
1) Mesh: to define the meshing between two gears.
2) Material: to define the material properties for a certain component (yield strength, ultimate tensile strength, etc.).
3) Force/Torque: to define any external forces or torques acting on a certain component.
4) Shaft Profile: to define the external profile of a shaft component.
5) Shaft Section: to define a section of a shaft whose profile has already been defined to analyze it.

The tool accounts for stress concentration factors, notch sensitivity factor, and fatigue limit correction factors to accurately perform the analysis.
The user must define the geometry of the shaft(s) properly so that the tool can account for all these factors accurately.

## Examples
**1. Normal Internal Load**

![Normal Load](Assets/normal_load.png)


**2. Bending Stress**

![Bending Stress](Assets/bending_stress.png)


**3. Haigh Diagram**

![Haigh Diagram](Assets/haigh_diagram.png)

## Installation
For installation, make sure the proper environment is activated and execute the following command in the command line prompt:
<pre> pip install pygritbx </pre>

## What's New?
- Bug fixes.
- Automated static and fatigue verification for shafts by implmenting the **".performStaticVerification()"** and **".performFatigueVerification()"** functions in the Shaft class.
- Automated gear tooth bending and gear tooth pitting analyses by implementing **".analyseGearToothBending()"** and **"analyseGearToothPitting()"** functions in the Gear class.
- Modified examples according to latest implementations.

## References
- Richard G. Budynas and J. Keith Nisbett, *Shigley's Mechanical Engineering Design*, McGraw-Hill, 2006.
- SKF Group, *SKF Rolling Bearings Catalogue*, SKF. [SKF Rolling Bearings](https://www.skf.com/group/products/rolling-bearings)