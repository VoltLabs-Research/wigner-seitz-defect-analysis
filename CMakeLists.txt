cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(wigner-seitz-defect-analysis VERSION 1.0.2 LANGUAGES CXX)

if(NOT TARGET coretoolkit::coretoolkit)
    find_package(coretoolkit REQUIRED)
endif()

volt_add_plugin(wigner-seitz-defect-analysis
    MAIN_SOURCE src/plugin.cpp
)
