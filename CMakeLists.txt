cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(wigner-seitz-defect-analysis VERSION 1.0.0 LANGUAGES CXX)

find_package(coretoolkit REQUIRED)

volt_add_plugin(wigner-seitz-defect-analysis
    MAIN_SOURCE src/plugin.cpp
)
