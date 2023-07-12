cmake_minimum_required(VERSION 3.25)
project(DG)

set(CMAKE_CXX_STANDARD 20)

add_executable(untitled main.cpp Gaussian.cpp DiscGalerkin.h Lagrangian.cpp GaussianLobatto.cpp)
