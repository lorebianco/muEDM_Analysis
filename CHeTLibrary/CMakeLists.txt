cmake_minimum_required(VERSION 3.16)
project(CHeTLibrary VERSION 1.0 LANGUAGES CXX)

# --- A. Smart Install Prefix ---
# If the user didn't specify a install path with -DCMAKE_INSTALL_PREFIX=...,
# default to a folder named "install" inside the project directory.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH "Default install path" FORCE)
endif()

# --- 0. Developer Experience (Editor Integration) ---
# Enable the generation of the compilation database
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# This creates a "task" that runs automatically during 'make'
# It copies the json file from build/ to the project root.
# 'copy_if_different' ensures we don't stress the disk if nothing changed.
add_custom_target(copy_compile_commands ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_BINARY_DIR}/compile_commands.json"
        "${CMAKE_SOURCE_DIR}/../compile_commands.json"
    COMMENT "Automatic: Copying compilation database to project root for Editor support"
    VERBATIM
)

# --- 1. Find ROOT Package ---
find_package(ROOT REQUIRED COMPONENTS Geom Gpad Graf3d RIO MathCore)
include(${ROOT_USE_FILE})

# --- 2. Define Source and Header Files ---
# Gather all .hh files from include/CHeT
file(GLOB_RECURSE HEADERS "include/CHeT/*.hh")
# Gather all .cc files from src/
file(GLOB_RECURSE SOURCES "src/*.cc")

# --- 3. Generate ROOT Dictionary ---
# This generates the G__CHeTLib.cxx file which allows ROOT to understand your custom classes
ROOT_GENERATE_DICTIONARY(G__CHeTLib
    ${HEADERS}
    LINKDEF include/CHeT/LinkDef.h
    OPTIONS -I${CMAKE_CURRENT_SOURCE_DIR}/include
)

# --- 4. Create Shared Library (.so) ---
# The library includes your source code plus the generated ROOT dictionary
add_library(CHeTLib SHARED ${SOURCES} G__CHeTLib.cxx)

# --- 5. Configure Include Directories ---
target_include_directories(CHeTLib PUBLIC
    # Interface for the current build environment
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    # Interface for installed/exported environment
    $<INSTALL_INTERFACE:include>
    # Include ROOT headers
    ${ROOT_INCLUDE_DIRS}
)

# --- 6. Link against ROOT Libraries ---
target_link_libraries(CHeTLib PUBLIC ${ROOT_LIBRARIES})

# --- 7. Install Rules ---
install(TARGETS CHeTLib DESTINATION lib)
install(DIRECTORY include/CHeT DESTINATION include)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/libCHeTLib_rdict.pcm"
    "${CMAKE_CURRENT_BINARY_DIR}/libCHeTLib.rootmap"
    DESTINATION lib
)
