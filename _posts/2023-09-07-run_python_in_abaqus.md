---
title: "How to generate FE model automatically by python in ABAQUS"
date: 2023-09-07
permalink: /posts/2023/09/python_abaqus/
tags:
  - python
  - ABAQUS
  - finite elements
---

In this post, we will explore how to use Python and Abaqus to create, simulate, and analyze 4D printed active composite structures. These structures have the unique ability to change their shape or properties in response to external stimuli, such as temperature changes.

# Generating and Analyzing 4D Printed Active Composite Structures Using Python

In this post, we will explore how to use Python and Abaqus to create, simulate, and analyze 4D printed active composite structures. These structures have the unique ability to change their shape or properties in response to external stimuli, such as temperature changes.

## Setting up the Environment

To begin, we need to set up our environment by importing the necessary modules and libraries:

```python
# -*- coding: mbcs -*-
# Import necessary modules and libraries
from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import random
```

## Defining Geometry
Next, we'll define the dimensions and material properties of our active composite structure. These parameters include:

* Dimensions of the geometry (length, width, cell size)
* Material properties for the 'active' and 'passive' materials (Young's Modulus, Poisson's Ratio, Coefficient of Thermal Expansion)

```python
# Define some variables for geometry

# Dimensions of the geometry
half_x = 80.0  # Half of the geometry's x-dimension
half_y = 10.0  # Half of the geometry's y-dimension
cell_size = 5.0  # Size of individual cells or elements in the mesh

# Material properties for the 'active' material
young_active = 2000000000.0  # Young's Modulus for the active material
poisson_active = 0.35  # Poisson's Ratio for the active material
expand_active = 0.001  # Coefficient of Thermal Expansion for the active material

# Material properties for the 'passive' material
young_passive = 1000000000.0  # Young's Modulus for the passive material
poisson_passive = 0.35  # Poisson's Ratio for the passive material
expand_passive = 0.0  # Coefficient of Thermal Expansion for the passive material

# Meshing parameters
seed_size = 0.5  # Size of the mesh seeds
seed_deviationFactor = 0.1  # Deviation factor for mesh seed size
seed_minSizeFactor = 0.1  # Minimum size factor for mesh seed size

# Temperature change for analysis (a tuple with one value)
temp_change = (100.0,)  # Temperature change applied in the analysis
```


Useful information
======

* Books
  * [Efficient Finite Element Modelling - Automated model generation & evaluation using Simulia Abaqus](http://Liuchao-JIN.github.io/files/research/abq_python.pdf)
* Videos
  * [How to create Python scripts automatically using ABAQUS CAE](https://www.youtube.com/watch?v=1zT2xGLv01E)
  * [Use coordinate instead of get Sequence from mask for material assign](https://www.youtube.com/watch?v=LVwKZnwh9W4)
  * [How to do ABAQUS Scripting; Simulating a Simple Disk Compression Test](https://www.youtube.com/watch?v=dKXgYoe3UPo)
  * [Learn ABAQUS Scripting; Export Results Automatically from ODB Files](https://www.youtube.com/watch?v=a5hT9oSPJKQ)
  * [How to READ and UNDERSTAND ABAQUS Files](https://www.youtube.com/watch?v=au7cQQU2rmE)
  * [How to get node labels labels in ABAQUS](https://www.youtube.com/watch?v=J3mzux9sQus)
* Websites
  * [Automatization of Abaqus FEA using python](https://balashov-artem.github.io/Portfolio/projects/automatization-of-abaqus-fea-using-python/)
  * [Using the Partition toolset in Abaqus](https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-c-parstart.htm)
  * [ABAQUS/CAE User's Manual](https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/usi/default.htm?startat=pt06ch44s04.html)
  * [Abaqus FEA Scripting with python](https://ifcuriousthenlearn.com/blog/2015/04/02/Abaqus-FEA-Scripting-with-python/)
