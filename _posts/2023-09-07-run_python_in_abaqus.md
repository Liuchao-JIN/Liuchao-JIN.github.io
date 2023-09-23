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


# Choosing Hardware for Abaqus: Factors to Consider


## Abaqus Runtime Speed Factors

When it comes to the performance of Abaqus, several factors come into play. We'll discuss two main aspects:

### 1. Factors Influencing Abaqus Runtime Speed

- **Hardware (70%)**
  - **CPU**
    - Number of Physical Cores: Abaqus supports multi-threading in Threads/MPI/Hybrid modes. Whether you're running it on a personal computer or a server cluster, more physical cores are better. However, there's a trade-off between cost and performance. The optimal core count for the best value is around 64. It's essential to note that we're referring to **physical cores** here, not logical threads, as explained in the next section. For personal use, it's recommended to have a minimum of 4 physical cores, with a preferable range of 8-32 cores. Avoid going beyond 64 cores as it's not cost-effective.
    - Memory Speed, Latency, and Total Bandwidth: More memory channels can improve processor access speed. Lower latency is better, and total bandwidth refers to the throughput of all memory channels used simultaneously. The memory performance for single-threaded tasks equals that of a single memory channel. When running multiple threads or processes, you can utilize multiple channels simultaneously. However, when multiple cores share a single memory channel, each core's memory performance decreases. Therefore, systems with the most memory channels and theoretically the highest bandwidth are optimal for memory-bound parallel applications like Abaqus.
    - Clock Frequency: A higher CPU clock frequency is preferable for Abaqus, regardless of the analysis type.
    - Mathematical Features (AVX2, AVX-512, FMA): The AVX instruction set is currently used only by the Abaqus/Standard solver. Having it can provide a slight speed boost (about 1-2% acceleration), but it's not essential.
    - Floating-Point Units (FPU): Data throughput depends on the clock frequency multiplied by the number of FPU units. For Abaqus, having additional FPUs is more valuable than increasing the clock frequency.
  - **RAM**
    - Match memory modules with CPU memory channels and total bandwidth. Choose higher-frequency RAM and fill the memory channels within your budget. The recommended core-to-memory ratio is 1:4 for the Explicit solver and 1:8 for the Standard solver. For example, if you have an 8-core CPU and primarily run Abaqus/Explicit, 32GB of RAM with 4 channels is a suitable choice.
  - **GPU**
    - For graphics cards, basic functionality is sufficient, as GPU acceleration is primarily limited to the direct sparse solver in Standard. Only a few GPU models support Abaqus acceleration, and they can be expensive. If budget allows and you frequently run Abaqus/Standard with the direct sparse solver, consider options like Tesla K40m, GP100/GV100, or AMD Radeon Pro V11.
  - **SSD/HDD**
    - It's advisable to install the operating system and software on a solid-state drive (SSD) and use a mechanical hard drive (HDD) for storing models and setting the working directory. Separating these components prevents interference, as Abaqus generates many temporary files (outside the working directory), which are automatically deleted after job completion. Having them on different drives avoids performance bottlenecks. If possible, consider using an SSD for the working directory as well. In tests, running a CEL airbag inflation model on an SSD-based working directory was approximately 12% faster than on an HDD.

  * **Note**: When considering hardware, also factor in the effects of hardware interconnectivity, which indirectly affect Abaqus' overall performance. Factors like cooling systems determine CPU throttling, and high-power CPUs benefit from water cooling. When building cluster servers, an increase in CPU count and MPI levels requires more MPI communication. Systems with higher bandwidth and lower latency interconnects often allow for greater performance scaling and higher CPU limits, especially when latency is significantly reduced, as is the case with InfiniBand interconnects compared to 10Gb Ethernet.

- **Operating System (3%)**
  - Differences in operating system types, versions, libraries, and configuration settings slightly affect Abaqus runtime speed. Linux commonly outperforms Windows in program scheduling, resulting in generally faster runtimes by up to 5%. However, for specific analysis types, it might be slightly slower. Therefore, the choice of the operating system is flexible and depends on personal preferences and needs.

- **Abaqus Task Properties (27%)**
  - Factors like model size, mesh characteristics, analysis procedures, solver types, precision, iteration or increment counts, and output settings all influence Abaqus calculation speed significantly. Many people overlook this aspect when selecting hardware for Abaqus, which can lead to slow performance on high-end systems. It's crucial to consider the following general guidelines when balancing hardware (mainly CPU and memory) performance:

    - Abaqus/Standard should have enough memory to run analyses, which is more critical than other factors. This is why the core-to-memory ratio is as high as 1:8, twice that of the Explicit solver.
    - For Abaqus/Explicit, the number of cores and memory access speed (bandwidth) are the primary factors.
    - For Abaqus/Standard direct solvers, massive models with continuous elements (e.g., dynamic systems or soil analysis) are computation-intensive, so higher clock frequencies are essential. Additionally, processor features like AVX2/AVX512 benefit from using Intel MKL for DGEMM calculations.
    - For Abaqus/Standard AMG iterative solvers or direct solvers with models not dominated by the solver (e.g., shell structures for aircraft fuselages), memory bandwidth is more critical.

## Abaqus: Intel or AMD CPU?

The answer to whether you should choose Intel or AMD processors for Abaqus depends largely on whether you primarily use Standard or Explicit solvers. In general:

- AMD offers a significant advantage in core count, providing better cost-effectiveness compared to Intel. If you prioritize core count and mainly perform calculations like those in the Explicit solver, AMD is the top choice, offering considerably higher performance than Intel processors in the same price range. For those on a budget, Ryzen processors are a good option. If you have a more substantial budget or your company is willing to invest, consider the Epyc second-generation (recommended 7F52/7H12) or third-generation (recommended 74F3/75F3/7T83) processors or the Threadripper series. These processors are designed specifically for CAE analysis, and AMD provides benchmarks for software like LS-DYNA and Abaqus, showcasing their excellent performance.
- Intel's advantage lies in its various instruction sets, which enhance the performance of Standard solvers and ensure perfect compatibility with software. AMD falls short in this regard, with lower versions of Abaqus experiencing compatibility issues on AMD platforms. For instance, even with Abaqus 2016.HF28 and later versions, there are still compatibility issues with running Standard and CFD co-simulation on AMD platforms. Since version 28 is likely the last maintenance patch, it means that the official support for this issue has ended. Therefore, if you are using a lower version of Abaqus (especially before 2018) and frequently use the Standard solver, be cautious when selecting AMD processors. In contrast, Intel processors do not face such problems. Recommended Intel processors include Core i7-12700K or i9-12900K, and using DDR5 memory with high clock frequencies can significantly improve Abaqus performance. For higher budgets, the cost-performance king is the 8375c, especially when paired with dual processors. This processor represents an outstanding choice for those who want optimal performance.

Finally, let's explain why we mentioned "physical cores" earlier. Intel's Hyper-Threading (HT) or AMD's Simultaneous Multithreading (SMT) technologies do not provide meaningful performance improvements for Abaqus parallel computing (in fact, they can even reduce performance for compute-intensive CAE analyses). These technologies only increase the number of logical threads and do not significantly boost performance. In recent years, Abaqus has recognized this fact, and starting from Abaqus 2022, it only recognizes "physical cores," regardless of whether HT or SMT is enabled in the BIOS. Users should specify the maximum number of "physical cores" to avoid errors. If you prefer to use HT/SMT and want to keep recognizing logical threads as in previous versions, you can add the following statement to the environment file (abaqus_v6.env):

```python
import os
os.environ['ABA_CPUS_LOGICAL'] = '1'
```
