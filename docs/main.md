# Galileo
*A light-weight and extensible C++ library for Guass-Legendre Pseudospectral Collocation of Switched Systems using Casadi and Pinocchio.*

<!-- <img src="https://i.imgur.com/VQJ3ZNe.png"/> -->

[![Documentation](https://img.shields.io/badge/docs-generate-brightgreen.svg)](https://github.com/echandler5956f/Galileo/tree/main/docs)
[![License BSD-3-Clause](https://img.shields.io/badge/license-MIT-blue.svg)](https://mit-license.org/)

Named after the famous scientist who posed one variation of the Brachistochrone problem, *Galileo* is an efficient optimal control framework that uses Gauss-Legendre Pseudospectral Collocation to solve the switched systems problem for legged robots.

Features:  

:heavy_check_mark: Intuitive and efficient formulation of variables, cost and constraints using [casadi].   

:heavy_check_mark: Solver interface enables using the high-performance solvers [Ipopt] and [Snopt].  

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

### From source with [CMake]

1. Install Galileo's mandatory dependencies:
   * [casadi]
   * [pinocchio]
   * [Eigen]
   * [Boost]

2. (optional) Install Galileo's optional dependencies
   * [Ipopt]

Once you have the required dependencies, you can install the library with
```bash
mkdir build && cd build && cmake .. && make -j4 && sudo make install
```

To uninstall the library, simply run
```bash
sudo make uninstall
```
from within the build directory.

## Run

We created a `test-installs` folder with some simple scripts to test your `casadi` + `pinocchio` install. Be sure to take a look if you are running into trouble. To test the actual library, please refer to the `examples` folder. If you have already built the repo from source, you can test it by running

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
in the main directory, and then
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

[casadi]: https://github.com/casadi/casadi
[pinocchio]: https://github.com/stack-of-tasks/pinocchio
[CMake]: https://cmake.org/cmake/help/v3.0
[ROS]: http://www.ros.org
[Ipopt]: https://projects.coin-or.org/Ipopt
[Snopt]: http://www.sbsi-sol-optimize.com/asp/sol_product_snopt.html
[rviz]: http://wiki.ros.org/rviz
[catkin]: http://wiki.ros.org/catkin
[catkin tools]: http://catkin-tools.readthedocs.org
[Eigen]: http://eigen.tuxfamily.org
[Boost]: https://www.boost.org