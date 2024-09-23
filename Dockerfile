FROM ros:noetic

ARG HSL_PATH=null

WORKDIR /home

ENV SHELL /bin/bash

# Set the default command to bash
CMD ["/bin/bash"]

# Set environment variables
ENV PATH="/usr/local/bin:$PATH"
ENV LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"
ENV CMAKE_PREFIX_PATH="/usr/local:$CMAKE_PREFIX_PATH"

# RUN apt-get update && apt-get upgrade -y
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
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
    libboost-all-dev \
    libassimp-dev \
    libclang-dev \
    llvm-dev \
    libblas-dev \
    libmetis-dev \
    libconsole-bridge-dev \
    liburdfdom-headers-dev \
    liboctomap-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN mkdir repos

# Eigen: pinocchio requirement
RUN cd repos && \
    curl -L https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz -o eigen-3.4.0.tar.gz && \
    tar -xzf eigen-3.4.0.tar.gz && rm eigen-3.4.0.tar.gz && \
    cd eigen-3.4.0 && \
    mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && make -j$(nproc) && make install && \
    cd /home && rm -rf repos/eigen-3.4.0

# HSL: casadi requirement
RUN cd repos && \
    git clone https://github.com/coin-or-tools/ThirdParty-HSL.git
ADD ${HSL_PATH} /home/repos/ThirdParty-HSL/coinhsl/
RUN cd repos/ThirdParty-HSL && \
    git checkout 4f8da75 && \
    rm -r .git && \
    ./configure && \
    make && \   
    make install && \
    cd /home && rm -rf repos/ThirdParty-HSL

# urdfdomm: pinocchio requirement
RUN cd repos && \
    git clone https://github.com/ros/urdfdom.git && \
    cd urdfdom && git checkout 3f6bf9a && rm -r .git && \
    mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && make -j$(nproc) && make install && \
    cd /home && rm -rf repos/urdfdom

# CppAD: CppADCodegen requirement
RUN cd repos && \
    git clone https://github.com/coin-or/CppAD.git && \
    cd CppAD && git checkout 470c769 && rm -r .git && \
    mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && make install && \
    cd /home && rm -rf repos/CppAD

# CppADCodegen: pinocchio requirement
RUN cd repos && \
    git clone https://github.com/joaoleal/CppADCodeGen.git && \
    cd CppADCodeGen && git checkout 656d23e && rm -r .git && \
    mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && make install && \
    cd /home && rm -rf repos/CppADCodeGen

# hpp-fcl: pinocchio requirement
RUN cd repos && \
    git clone https://github.com/humanoid-path-planner/hpp-fcl.git && \
    cd hpp-fcl && git checkout 6b9f9c8 && rm -r .git && \
    mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_INTERFACE=OFF -DHPP_FCL_HAS_OCTOMAP=ON && \
    make -j$(nproc) install && \
    cd /home && rm -rf repos/hpp-fcl

# casadi: galileo, pinocchio requirement
RUN cd repos && \
    git clone https://github.com/casadi/casadi.git && \
    cd casadi && git checkout 81bbcd3 && rm -r .git && \
    mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DWITH_BUILD_IPOPT=ON -DWITH_IPOPT=ON -DWITH_HSL=ON -DWITH_OPENMP=ON -DWITH_MUMPS=ON -DWITH_BUILD_MUMPS=ON -DWITH_SNOPT=OFF -DWITH_CLANG=ON -DWITH_THREAD=ON -DWITH_OSQP=ON -DWITH_BUILD_OSQP=ON -DWITH_QPOASES=ON -DWITH_LAPACK=ON -DWITH_BUILD_LAPACK=ON -DWITH_BLOCKSQP=ON -DWITH_PYTHON=OFF -DWITH_PYTHON3=OFF -DWITH_BUILD_METIS=ON .. && \
    make -j$(nproc) && make install && \
    cd /home && rm -rf repos/casadi

# pinocchio: galileo requirement
RUN cd repos && \
    git clone --recursive https://github.com/stack-of-tasks/pinocchio && \
    cd pinocchio && \ 
    git checkout 0b594a0 && \
    rm -r .git && \
    mkdir build && \ 
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_BENCHMARK=ON -DBUILD_UTILS=ON -DBUILD_PYTHON_INTERFACE=OFF -DGENERATE_PYTHON_STUBS=OFF -DBUILD_WITH_URDF_SUPPORT=ON -DBUILD_WITH_COLLISION_SUPPORT=ON -DBUILD_WITH_AUTODIFF_SUPPORT=ON  -DBUILD_WITH_CASADI_SUPPORT=ON -DBUILD_WITH_CODEGEN_SUPPORT=ON -DBUILD_WITH_OPENMP_SUPPORT=ON .. && \
    make -j$(nproc) && \
    make install && \
    cd /home && rm -rf repos/pinocchio

# Galileo
# RUN cd repos && \
#     git clone https://github.com/KaiNakamura/Galileo.git --depth=1 && \
#     cd Galileo && \ 
#     mkdir build && \ 
#     cd build && \
#     cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_WITH_OPENMP=ON .. && \
#     make && \
#     make install
