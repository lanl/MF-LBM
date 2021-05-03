<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Thanks again! Now go create something AMAZING! :D
-->

# MF-LBM: A Portable, Scalable and High-performance Lattice Boltzmann Code for Flow in Porous Media


<img src="images/core-dns.png" alt="drawing" width="900"/>

<br/>
<br/>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-code">About The Code</a></li>
        <ul>
          <li><a href="#features">Features</a></li>
        </ul>
        <ul>
          <li><a href="#components">Components</a></li>
        </ul>
        <ul>
          <li><a href="#technical-details-of-the-main-simulation-code">Technical details of the main simulation code</a></li>
        </ul>
        <ul>
          <li><a href="#citing-mf-lbm">Citing MF-LBM</a></li>
        </ul>
    <li><a href="#build-instructions">Build Instructions</a></li>
        <ul>
          <li><a href="#prerequisites">Prerequisites</a></li>
        </ul>
        <ul>
          <li><a href="#installation">Installation</a></li>
        </ul>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#important-notes">Important Notes</a></li> 
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#references">References</a></li>
  </ol>
</details>

<br/>
<!-- ABOUT THE CODE -->

## About The Code

MF-LBM (Chen et al., 2018, 2019) is a high-performance lattice Boltzmann (LB) code for direct numerical simulation (DNS) of flow in porous media, primarily developed by Dr. Yu Chen (LANL), under the supervision of Prof. Albert Valocchi (UIUC), Dr. Qinjun Kang (LANL) and Dr. Hari Viswananthan (LANL). 'MF' refers to microfluidics or 'Magic Find'. The code was first developed at University of Illinois at Urbana-Champaign based on a mainstream LB color-gradient multiphase model and further improved at Los Alamos National Laboratory by implementing the Continuum-Surface-Force and geometrical wetting models to reduce spurious currents so that the inertial effects in scCO<sub>2</sub> and brine displacement in porous media can be accounted for (Chen et al., 2019). 

### Features
* exploring multiple levels of parallelism
* extensively optimized for vectorization
* directive-based parallel programing model supporting CPU, GPU, MIC and ARM
* advanced LB multiphase model (CSF model + geometrical wetting model) ensuring relatively small spurious currents
* overlapped communication and computation
* pre-processing and post-processing code included
  
### Components
* Pre-processing code: 
  1) converting text images from the scans to a single 3D wall array stored in binary formate. 
  2) modifying the 3D wall array, i.e., adding inlet/outlet buffer layers 
  3) obtaining orientations of the solid surface stored in binary formate.
* Main simulation code:
  1) multiphase flow simulation code
  2) single phase flow simulation code (comming soon)
* Post-processing code:
  1) converting distributed data files to single vtk file
  2) reading in distributed data files for further analysis
### Technical details of the main simulation code
Modern manycore processors/coprocessors, such as GPUs and Intel Xeon Phi processors, are developing rapidly and greatly boost computing power. These processors not only provide much higher FLOPS (floating-point operations per second) but also much higher memory bandwidth compared with traditional CPU. One of the most attractive features of the lattice Boltzmann method (LBM) is that it is explicit in time, has nearest neighbor communication, and the computational effort is in the collision step, which is localized at a grid node. For these reasons, the LBM is well suited for manycore processors which require a higher degree of explicit parallelism. The data movement in the LBM is much more intensive than for traditional CFD considering that the D3Q19 lattice model has 19 lattice velocities. Given the current state of computational hardware, in particular the relative speed and capacity of processors and memory, the LBM is a memory-bandwidth-bound numerical method. The high memory bandwidth provided by GPUs or Intel Xeon Phi processors greatly benefits the LBM.

The code is written on Fortran 90 and employs MPI-OpenACC/OpenMP hybrid programing model. The main reason that we chose OpenACC/OpenMP (directive-based parallel programming models) over CUDA is that we want to keep the code portable across different computing platforms so that we are not limited by the NVIDIA GPU solution. As GPU and Intel Xeon Phi processor (and even latest CPU from Intel with AVX512 instructions) rely heavily on SIMT/SIMD, the optimization strategy for these manycore processors/coprocessors are similar, which enables us to achieve reasonable performance across different platforms:
* The AA pattern streaming method is employed to significantly reduce memory access and memory consumption.
* The structure of arrays (SoA) data layout is used to achieve coalesced memory access and maximize vectorization.
* Communication and computation is overlapped to achieve good parallel efficiency, particularly for heterogenous computing platforms.

### Citing MF-LBM
Chen, Y., Valocchi, A., Kang, Q., & Viswanathan, H. S. (2019). Inertial effects during the process of supercritical CO2 displacing brine in a sandstone: Lattice Boltzmann simulations based on the continuum‐surface‐force and geometrical wetting models. Water Resources Research, 55, 11144– 11165. https://doi.org/10.1029/2019WR025746

<br/>

<!-- Build Instructions -->
## Build Instructions

### Prerequisites
* A Fortran compiler 
* A MPI implementation
* CUDA toolkit (for NVIDIA GPU platform)
* PGI Fortran compiler (for NVIDIA GPU platform)
* Intel Fortran compiler (for Intel Xeon Phi)
* Make

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/your_username_/Project-Name.git
   ```
2. Make
   1. CPU version
        ```sh
        cd mf-lbm/3d-multiphase
        # Make necessary changes to makefile.
        # Choose cpu as the architecture option.
        # Choose compiler.
        # OpenMP is recommended for CPU or ARM version.
        # See instruction in makefile for more information.  
        your-preferred-editor makefile
        make
        ```
    1. GPU version
        ```sh
        cd mf-lbm/3d-multiphase
        # Make necessary changes to makefile.
        # Choose gpu as the architecture option.
        # Choose compiler.
        # OpenMP must be disabled for GPU version
        # See instruction in makefile for more information.  
        your-preferred-editor makefile
        # Make sure to enable OpenACC in preprocessor  
        your-preferred-editor mf-lbm/3d-multiphase/0.src/preprocessor.h  
 
        make
        ```
    2. MIC (Intel Xeon Phi) version
        ```sh
        cd mf-lbm/3d-multiphase
        # Make necessary changes to makefile.
        # Choose mic as the architecture option.
        # Choose compiler.
        # OpenMP and AVX512 must be enabled for MIC version
        # See instruction in makefile for more information.  
        your-preferred-editor makefile
        your-preferred-editor makefile     
        make
        ```
3. Configure run scripts (for CPU platform)
      ```sh
      cd working_directory
      cp path-to-repo/3d-multiphase/run_template/template-config_sim.sh ./config_sim.sh
      # Make necessary changes to config_sim.sh (i.e., input parameters, paths and run command). 
      # See instruction in config_sim.sh for more information.
      your-preferred-editor config_sim.sh
      ./config_sim.sh
      # The run-program script, irun.sh, will be generated. 
      
      # If OpenMP is enabled (recommended for CPU, MIC, and ARM platform); then run following
      # export OMP_NUM_THREADS=n
      # here n is recommended to be the core count of the NUMA of the CPU.

      # MPI process number is recommended to be the total available number of NUMA.

      # In the GPU version, OpenMP is disabled while OpenAcc is enabled.
      ```
4. Run the program
   ```sh
   ./irun.sh new
   ```
<br/>

<!-- USAGE EXAMPLES -->
## Usage
### The main simulation code:

* A drop attached to a flat wall
   ```sh
   cd working_directory
   cp path-to-repo/3d-multiphase/test_suites/3D_simulation/1.drop_attached_wall/config.sh ./
   # necessary modifications might needed according to your computing environment
   your-preferred-editor config_sim.sh
   ./config_sim.sh    
   ./irun.sh new
   ```
* Nonwetting fluid displacing wetting fluid in a square duct
   ```sh
   cd working_directory
   cp path-to-repo/3d-multiphase/test_suites/3D_simulation/2.drainage/config.sh ./
   # necessary modifications might needed according to your computing environment
   your-preferred-editor config.sh
   ./config_sim.sh    
   ./irun.sh new
   ```
* Nonwetting fluid displacing wetting fluid in a square duct with hard-coded geometry
   ```sh
   cd working_directory
   cp path-to-repo/3d-multiphase/test_suites/3D_simulation/3.drainage_external_geometry/config.sh ./
   # necessary modifications might needed according to your computing environment
   your-preferred-editor config_sim.sh
   ./config_sim.sh    
   ./irun.sh new
   ```
* Wetting fluid displacing nonwetting fluid in a square duct with external geometry file
   ```sh
   cd working_directory
   cp path-to-repo/3d-multiphase/test_suites/3D_simulation/3.drainage_external_geometry/config.sh ./
   # necessary modifications might needed according to your computing environment
   # specify the geometry file path on config.sh
   your-preferred-editor config.sh
   ./config_sim.sh    
   ./irun.sh new
   ```
* Body force driven fractional flow in a square duct with external geometry file and pre-processed solid-boudnary-info file
   ```sh
   cd working_directory
   cp path-to-repo/3d-multiphase/test_suites/3D_simulation/3.drainage_external_geometry/config.sh ./
   # necessary modifications might needed according to your computing environment
    # specify the geometry file and solid-boundary-info file paths on config.sh
   your-preferred-editor config.sh
   ./config_sim.sh    
   ./irun.sh new
   ```

### Pre-processing code:

### Post-processing code:

<br/>

## Important Notes
* Geometry preprocessing
* Contact angle
* Run command
* Domain decomposition
* Number of threads in OpenMP
* 
* GCC10 compiler issue: 
  Building the code with GCC10 may show error messages like 
  >Type mismatch between actual argument at (1) and actual argument at (2)

  This is a known issue with GCC10. Use GCC10 new option *-fallow-argument-mismatch* to turn these errors to warnings.

<br/>

<!-- LICENSE -->
## License
Distributed under the BSD-3 License. See [LICENSE] for more information.

© 2021. Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

<br/>

<!-- CONTACT -->
## Contact

Dr. Yu Chen - yu_chen_007@outlook.com

Dr. Qinjun Kang - qkang@lanl.gov

<br/>

## References

1. Chen, Yu, Albert J. Valocchi, Qinjun Kang, and Hari S. Viswanathan. "Inertial effects during the process of supercritical CO2 displacing brine in a sandstone: Lattice Boltzmann simulations based on the continuum‐surface‐force and geometrical wetting models." Water Resources Research 55, no. 12 (2019): 11144-11165.
2. Chen, Yu, Yaofa Li, Albert J. Valocchi, and Kenneth T. Christensen. "Lattice Boltzmann simulations of liquid CO2 displacing water in a 2D heterogeneous micromodel at reservoir pressure conditions." Journal of contaminant hydrology 212 (2018): 14-27.



<!-- MARKDOWN LINKS & IMAGES -->

[LICENSE]: https://github.com/ychen-hpc/mf-lbm-dev/blob/master/LICENSE

