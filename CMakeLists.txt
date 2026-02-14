cmake_minimum_required(VERSION 3.15)
project(ATPGeomodeling LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Library from src
add_library(atpgeomodeling STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/CATPGeoModel.cpp
)
target_include_directories(atpgeomodeling
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Eigen3 (header-only)
find_package(Eigen3 QUIET CONFIG)
if(Eigen3_FOUND)
    target_link_libraries(atpgeomodeling PUBLIC Eigen3::Eigen)
else()
    find_package(Eigen3 QUIET)
    if(Eigen3_FOUND)
        target_link_libraries(atpgeomodeling PUBLIC Eigen3::Eigen)
    else()
        message(WARNING "Eigen3 not found. Set EIGEN3_INCLUDE_DIR to your Eigen headers if needed.")
        if(DEFINED EIGEN3_INCLUDE_DIR)
            target_include_directories(atpgeomodeling PUBLIC ${EIGEN3_INCLUDE_DIR})
        endif()
    endif()
endif()

# vcpkg convenience: auto-detect VTK_DIR from VCPKG_ROOT if possible
if(NOT DEFINED VTK_DIR AND DEFINED ENV{VCPKG_ROOT})
    set(_vcpkg_root "$ENV{VCPKG_ROOT}")
    if(NOT DEFINED VCPKG_TARGET_TRIPLET)
        if(CMAKE_SIZEOF_VOID_P EQUAL 8)
            set(_vcpkg_triplet "x64-windows")
        else()
            set(_vcpkg_triplet "x86-windows")
        endif()
    else()
        set(_vcpkg_triplet "${VCPKG_TARGET_TRIPLET}")
    endif()
    set(_vcpkg_prefix "${_vcpkg_root}/installed/${_vcpkg_triplet}")
    if(EXISTS "${_vcpkg_prefix}/share/vtk/vtk-config.cmake")
        set(VTK_DIR "${_vcpkg_prefix}/share/vtk")
        list(PREPEND CMAKE_PREFIX_PATH "${_vcpkg_prefix}")
    endif()
endif()

# VTK
find_package(VTK CONFIG REQUIRED)

# Link VTK to the library so its usage requirements propagate
if(DEFINED VTK_LIBRARIES)
    target_link_libraries(atpgeomodeling PUBLIC ${VTK_LIBRARIES})
endif()

# Executable from example
add_executable(ATPExtrudeWithConstrain
    ${CMAKE_CURRENT_SOURCE_DIR}/example/ATPExtrudeWithConstrain.cpp
)

target_link_libraries(ATPExtrudeWithConstrain
    PRIVATE
        atpgeomodeling
)

set(_example_data
    ${CMAKE_CURRENT_SOURCE_DIR}/example/surf_frac_3.neu
    ${CMAKE_CURRENT_SOURCE_DIR}/example/surf2.stl
    ${CMAKE_CURRENT_SOURCE_DIR}/example/surf3.stl
)
foreach(_f IN LISTS _example_data)
    add_custom_command(TARGET ATPExtrudeWithConstrain POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${_f}" "$<TARGET_FILE_DIR:ATPExtrudeWithConstrain>"
    )
endforeach()

# If modern VTK targets are available, also link them explicitly and autoinit
if(TARGET VTK::CommonCore)
    # Collect commonly used VTK modules; add only those that exist
    set(_vtk_candidates
        VTK::CommonCore
        VTK::CommonDataModel
        VTK::FiltersCore
        VTK::FiltersGeneral
        VTK::FiltersSources
        VTK::IOXML
        VTK::IOLegacy
        VTK::IOGeometry
        VTK::InteractionStyle
        VTK::IOGAMBIT
    )
    set(_vtk_targets "")
    foreach(_m IN LISTS _vtk_candidates)
        if(TARGET ${_m})
            list(APPEND _vtk_targets ${_m})
        endif()
    endforeach()
    if(_vtk_targets)
        target_link_libraries(atpgeomodeling PUBLIC ${_vtk_targets})
        target_link_libraries(ATPExtrudeWithConstrain PRIVATE ${_vtk_targets})
    endif()
endif()

# VTK module auto-init (needed for object factory in VTK >= 8.90)
if(COMMAND vtk_module_autoinit)
    vtk_module_autoinit(
        TARGETS atpgeomodeling ATPExtrudeWithConstrain
        MODULES ${VTK_LIBRARIES}
    )
endif()

# Organize output folders (optional, useful on multi-config generators)
foreach(_t atpgeomodeling ATPExtrudeWithConstrain)
    set_target_properties(${_t} PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
    foreach(_cfg Debug Release RelWithDebInfo MinSizeRel)
        string(TOUPPER "${_cfg}" _ucfg)
        set_target_properties(${_t} PROPERTIES
            ARCHIVE_OUTPUT_DIRECTORY_${_ucfg} "${CMAKE_BINARY_DIR}/lib"
            LIBRARY_OUTPUT_DIRECTORY_${_ucfg} "${CMAKE_BINARY_DIR}/lib"
            RUNTIME_OUTPUT_DIRECTORY_${_ucfg} "${CMAKE_BINARY_DIR}/bin"
        )
    endforeach()
endforeach()
