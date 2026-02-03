cmake_minimum_required(VERSION 3.5)

project(smallpad LANGUAGES CXX)
enable_language(Fortran)


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(RES_FILES "")

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package (Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/lapack)
add_library(dcl libamos.h libamos.f libeispack.cpp libeispack.h libdcl.cpp libdcl.h)
target_link_libraries(dcl lapack blas gfortran m Eigen3::Eigen)
find_package(Qt5PrintSupport REQUIRED)

if(WIN32)
    find_package(QT NAMES Qt6 Qt5 COMPONENTS Widgets REQUIRED)
    find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets REQUIRED)
else()
    find_package(QT NAMES Qt6 Qt5 COMPONENTS Widgets X11Extras REQUIRED)
    find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets X11Extras REQUIRED)
endif()

if(WIN32)
 set(RES_FILES "resource.rc")
 set(CMAKE_RC_COMPILER_INIT windres)
 ENABLE_LANGUAGE(RC)
 SET(CMAKE_RC_COMPILE_OBJECT "<CMAKE_RC_COMPILER> <FLAGS> -O coff <DEFINES> -i <SOURCE> -o <OBJECT>")
    set(PROJECT_SOURCES
        main.qrc
        tinyexpr/tinyexpr.cpp
        tinyexpr/tinyexpr.h
        main.cpp
        qcustomplot.cpp
        qcustomplot.h
        pad.cpp
        pad.h
        pad.ui
        about.ui
        settings.ui
        integrate.ui
        plot.ui
        solve.ui
        check.ui
        )
else()
    set(PROJECT_SOURCES
        main.qrc
        tinyexpr/tinyexpr.cpp
        tinyexpr/tinyexpr.h
        main.cpp
        qcustomplot.cpp
        qcustomplot.h
        pad.cpp
        nativeeventfilter.h
        nativeeventfilter.cpp
        pad.h
        pad.ui
        about.ui
        settings.ui
        integrate.ui
        plot.ui
        solve.ui
        check.ui
        )
endif()
if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(smallpad
        ${PROJECT_SOURCES}
        ${RES_FILES}
        ${CMAKE_CURRENT_LIST_DIR}/main.qrc
        )
else()
    if(ANDROID)
        add_library(smallpad SHARED
            ${PROJECT_SOURCES}
            )
    else()
        add_executable(smallpad
            ${PROJECT_SOURCES}
            ${RES_FILES}
            ${CMAKE_CURRENT_LIST_DIR}/main.qrc
            )
    endif()
endif()
if(WIN32)
    set_property(TARGET smallpad PROPERTY WIN32_EXECUTABLE true)
    set(app_icon_resource_windows "${CMAKE_CURRENT_LIST_DIR}/resource.rc")
    target_link_libraries(smallpad PRIVATE Qt${QT_VERSION_MAJOR}::Widgets Qt5::PrintSupport user32 dcl)
else()
    target_link_libraries(smallpad PRIVATE Qt${QT_VERSION_MAJOR}::Widgets X11 Xtst Qt5::PrintSupport Qt5::X11Extras dcl xdo)
endif()
