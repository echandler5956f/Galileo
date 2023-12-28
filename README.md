# Galileo
*A light-weight and extensible C++ library for Guass-Legendre Pseudospectral Collocation of Switched Systems using CasADi and Pinocchio.*

<img src="https://i.imgur.com/VQJ3ZNe.png"/>

[![Documentation](https://img.shields.io/badge/docs-generate-brightgreen.svg)](https://github.com/echandler5956f/Galileo/tree/main/docs)
[![License BSD-3-Clause](https://img.shields.io/badge/license-MIT-blue.svg)](https://mit-license.org/)

Named after the famous scientist who posed one variation of the Brachistochrone problem, *Galileo* is an efficient optimal control framework that uses Gauss-Legendre Pseudospectral Collocation to solve the switched systems problem for legged robots.

Features:  

:heavy_check_mark: Intuitive and efficient formulation of variables, cost and constraints using [CasADi].   

:heavy_check_mark: Solver interface enables using the high-performance solvers [Ipopt] and [SNOPT].  

:heavy_check_mark: [pinocchio] makes custom robot integration as simple as switching the URDF.

:heavy_check_mark: [ROS]/[catkin] integration (optional).

:heavy_check_mark: Light-weight framework makes it easy to use and extend.

<br>

<p align="center">
  <a href="#install">Install</a> •
  <a href="#run">Run</a> •
  <a href="#develop">Develop</a> •
  <a href="#contribute">Contribute</a> •
  <a href="#publications">Publications</a> •
  <a href="#credits">Credits</a>
</p>
<br/>

## Install

The following Linux installation instructions are provided for your convenience:

### From source with [CMake]

1. Install Galileo's mandatory dependencies:
   * [CasADi]
   * [pinocchio]
   * [Eigen]
   * [Boost]

2. Install Galileo's optional dependencies
   * [SNOPT]
   * [HSL]
   * [OpenMP]
   * [gnuplot]

### CasADi Source Install
First, gather the dependencies listed here:
https://github.com/casadi/casadi/wiki/InstallationLinux

Then, clone the repo and make a build folder with
```bash
git clone https://github.com/casadi/casadi.git && cd casadi && mkdir build && cd build
```

Now, run
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DWITH_BUILD_IPOPT=ON -DWITH_IPOPT=ON -DWITH_MUMPS=ON -DWITH_BUILD_MUMPS=ON  ..
```
If you have an [HSL] license (highly recommended, as these solvers tend to speed up convergence by ~2x for our problems, and academics can acquire one for free!), you should manually build the [ThirdParty-HSL] interface and then add the following CasADi flag:

```bash
 -DWITH_HSL=ON
 ```

 If you want to compile [CasADi] with [OpenMP] (recomended), you can add
 ```bash
 -DWITH_OPENMP=ON
 ```

 Similarly, you can add the optional [SNOPT] interface with
 ```bash
 -DWITH_SNOPT=ON
 ```

 Now, build [CasADi] from source:
 ```bash
 make
 sudo make install
 ```

### Pinocchio Source Install

First, gather the dependencies listed here:
https://stack-of-tasks.github.io/pinocchio/download.html

Clone the repo and its submodules:
```bash
git clone --recursive https://github.com/stack-of-tasks/pinocchio
```

Make a build directory
```bash
cd pinocchio && mkdir build && cd build
```

Run cmake with [CasADi] support 
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_WITH_CASADI_SUPPORT=ON ..
```

and then build and install with
```bash
make -j4
sudo make install
```

### Gnuplot Optional Install
Gnuplot is our plotting library of choice. Installation is very straightforward from the package manager with:
```bash
sudo apt-get update
sudo apt-get install gnuplot
```

### Galileo Source Install

Once you have the required dependencies for Galileo, clone the repo
```bash
git clone https://github.com/echandler5956f/Galileo.git
```

create a build directory
```bash
cd Galileo && mkdir build && cd build
```

and finally, build and install Galileo with
```bash
 cmake -DCMAKE_BUILD_TYPE=Release .. && make && sudo make install
```

Optionally, you can add the `-DBUILD_WITH_OPENMP=ON` flag to enable using [OpenMP] for parallel evaluation of the constraint maps (highly recommended). Note that you must have enabled the OpenMP interface when installing [CasADi] for this to work.

To uninstall the library, simply run
```bash
sudo make uninstall
```
from within the build directory.

## Run

We created a `test-installs` folder with some simple scripts to test your [CasADi] + [pinocchio] install. Be sure to take a look if you are running into trouble. To test the actual library, please refer to the `examples` folder. If you have already built the repo from source, you can test it by running

```bash
build/examples/simple_test
```
or
```bash
build/examples/model_building_test
```
or
```bash
build/examples/traj_test
```
in the main repo directory.

## Documentation

### Doxygen
Run 
```bash
doxygen
```
in the main directory, and then open the `index.html` file in the `docs` folder. If you are using WSL like me, you can run
```bash
cd docs/html && explorer.exe index.html
```
to view the doxygen output.

## Develop

Right now, our main development goal is to get the framework to solve general trajectory optimization problems with fixed contact sequences in an MPC context. This is essentially the same problem that OCS2 and Crocoddyl solve. Our hope is that using Pseudospectral Collocation will yield faster convergence, and that the Casadi + Pinocchio pipeline will result in smaller, more digestible code, which is easier to expand upon than the bulky frameworks mentioned prior.

## Contribute

We are a small team of WPI students with no lab and no funding. We are open to contributions from the community! Submit a pull request or post an issue if you have any suggestions. This framework is still in its most nascent phase, and will take time to mature. Be patient and hopefully something good will come from this work.

## Publications

Coming soon (ICRA 2025?).

## Credits

### Written by 

- Ethan Chandler
- Akshay Jaitly

### With contributions from

- Hushmand Esmaeili
- Ibrahim Salman Al-Tameemi
- Duc Doan
- Zhun Cheng
<!-- - Yifu Yuan -->
<!-- - Lehong Wang -->
<!-- - Puen Xu -->
<!-- - Tao Zou -->
<!-- - Nikhil Gangaram -->
<!-- - Dheeraj Bhogisetty -->
<!-- - Nhi Nguyen -->

[CasADi]: https://github.com/casadi/casadi
[pinocchio]: https://github.com/stack-of-tasks/pinocchio
[CMake]: https://cmake.org/cmake/help/v3.0
[ROS]: http://www.ros.org
[Ipopt]: https://projects.coin-or.org/Ipopt
[SNOPT]: http://www.sbsi-sol-optimize.com/asp/sol_product_snopt.html
[rviz]: http://wiki.ros.org/rviz
[catkin]: http://wiki.ros.org/catkin
[catkin tools]: http://catkin-tools.readthedocs.org
[Eigen]: http://eigen.tuxfamily.org
[Boost]: https://www.boost.org
[gnuplot]: https://sourceforge.net/p/gnuplot/gnuplot-main/ci/master/tree/
[OpenMP]: https://www.openmp.org/
[HSL]: https://www.hsl.rl.ac.uk/
[ThirdParty-HSL]: https://github.com/coin-or-tools/ThirdParty-HSL