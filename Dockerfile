FROM ros:noetic


WORKDIR /home

ENV SHELL /bin/bash

# Set the default command to bash
CMD ["/bin/bash"]

# RUN apt-get update && apt-get upgrade -y
RUN apt-get update
RUN apt install -y \
    gcc \
    g++ \
    gfortran \
    git \
    cmake \
    wget \
    make \
    liblapack-dev \
    pkg-config \
    gnuplot \
    libeigen3-dev \
    libboost-all-dev \
    libassimp-dev \
    libclang-dev \
    llvm-dev \
    libblas-dev \
    libmetis-dev \
    libconsole-bridge-dev \
    liburdfdom-headers-dev \
    --install-recommends

ARG HSL_PATH=null

RUN mkdir repos
# HSL
RUN cd repos && \
    git clone https://github.com/coin-or-tools/ThirdParty-HSL.git
ADD ${HSL_PATH} /home/repos/ThirdParty-HSL/coinhsl/
RUN cd repos/ThirdParty-HSL && \
    git checkout 4f8da75 && \
    rm -r .git && \
    ./configure && \
    make && \   
    make install

# casadi
RUN cd repos && \
    git clone https://github.com/casadi/casadi.git && \ 
    cd casadi && \ 
    git checkout 81bbcd3 && \
    rm -r .git && \
    mkdir build && \ 
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DWITH_BUILD_IPOPT=ON -DWITH_IPOPT=ON -DWITH_HSL=ON -DWITH_OPENMP=ON -DWITH_MUMPS=ON -DWITH_BUILD_MUMPS=ON -DWITH_SNOPT=OFF -DWITH_CLANG=ON -DWITH_THREAD=ON -DWITH_OSQP=ON -DWITH_BUILD_OSQP=ON -DWITH_QPOASES=ON -DWITH_LAPACK=ON -DWITH_BUILD_LAPACK=ON -DWITH_BLOCKSQP=ON -DWITH_PYTHON=OFF -DWITH_PYTHON3=OFF -DWITH_BUILD_METIS=ON .. && \
    # cmake -DCMAKE_BUILD_TYPE=Release -DWITH_BUILD_IPOPT=ON -DWITH_IPOPT=ON -DWITH_HSL=ON -DWITH_OPENMP=ON -DWITH_MUMPS=ON -DWITH_BUILD_MUMPS=ON -DWITH_SNOPT=ON -DWITH_CLANG=ON -DWITH_THREAD=ON -DWITH_OSQP=ON -DWITH_BUILD_OSQP=ON -DWITH_QPOASES=ON -DWITH_LAPACK=ON -DWITH_BUILD_LAPACK=ON -DWITH_BLOCKSQP=ON -DWITH_PYTHON=ON -DWITH_PYTHON3=ON .. && \
    make && \
    make install

# urdfdomm
RUN cd repos && \
    git clone https://github.com/ros/urdfdom.git && \
    cd urdfdom && \ 
    git checkout 3f6bf9a && \
    rm -r .git && \
    mkdir build && \ 
    cd build && \
    cmake ../ && \
    make -j8 && \
    make install

# pinocchio
# ENV urdfdom_headers_DIR=“usr/include”
RUN cd repos && \
    git clone --recursive https://github.com/stack-of-tasks/pinocchio && \
    cd pinocchio && \ 
    git checkout 0b594a0 && \
    rm -r .git && \
    mkdir build && \ 
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_BENCHMARK=ON -DBUILD_UTILS=ON -DBUILD_PYTHON_INTERFACE=OFF -DGENERATE_PYTHON_STUBS=OFF -DBUILD_WITH_URDF_SUPPORT=ON -DBUILD_WITH_COLLISION_SUPPORT=ON -DBUILD_WITH_AUTODIFF_SUPPORT=ON  -DBUILD_WITH_CASADI_SUPPORT=ON -DBUILD_WITH_CODEGEN_SUPPORT=ON -DBUILD_WITH_OPENMP_SUPPORT=ON .. && \
    make -j4 && \
    make install

# CppAD
RUN cd repos && \
    git clone https://github.com/coin-or/CppAD.git && \
    cd CppADCodeGen && \ 
    git checkout 470c769 && \
    rm -r .git && \
    mkdir build && \ 
    cd build && \
    cmake .. && \
    make install

# CppADCodegen 
RUN cd repos && \
    git clone https://github.com/joaoleal/CppADCodeGen.git && \
    cd CppADCodeGen && \ 
    git checkout 656d23e && \
    rm -r .git && \
    mkdir build && \ 
    cd build && \
    cmake ../CppADCodeGen && \
    make install

# Galileo
RUN cd repos && \
    git clone https://github.com/echandler5956f/Galileo.git --depth=1 && \
    cd Galileo && \ 
    mkdir build && \ 
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_WITH_OPENMP=ON .. && \
    make && \
    make install
