cmake_minimum_required(VERSION 3.14) # Versione minima per FetchContent moderno

project(myfitter)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- 1. FETCHCONTENT: FTXUI ---
include(FetchContent)
FetchContent_Declare(ftxui
  GIT_REPOSITORY https://github.com/ArthurSonzogni/FTXUI
  GIT_TAG v5.0.0
)
FetchContent_MakeAvailable(ftxui)

# --- 2. FIND PACKAGES ---
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist MathCore Geom Eve)
include_directories(${ROOT_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})
include_directories(${Geant4_INCLUDE_DIRS})

# --- 3. CHET LIBRARY SETUP ---
set(CHeTLibrary_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../CHeTLibrary/install)
include_directories(${CHeTLibrary_DIR}/include)
link_directories(${CHeTLibrary_DIR}/lib)

# --- 4. GENFIT SETUP ---
if(DEFINED ENV{GENFIT})
    set(GENFIT_DIR $ENV{GENFIT})
    message(STATUS "Using GenFit from ${GENFIT_DIR}")
else()
    message(FATAL_ERROR "GENFIT environment variable not set.")
endif()

include_directories(
    ${GENFIT_DIR}/include
    ${GENFIT_DIR}/core/include
    ${GENFIT_DIR}/fitters/include
    ${GENFIT_DIR}/measurements/include
    ${GENFIT_DIR}/fields/include
    ${GENFIT_DIR}/trackReps/include
    ${GENFIT_DIR}/utilities/include
    ${GENFIT_DIR}/eventDisplay/include
    ${GENFIT_DIR}/GBL/include
)
link_directories(${GENFIT_DIR}/lib)

# --- 5. PROJECT INCLUDES ---
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/external)

# --- 6. SOURCES ---
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# --- 7. EXECUTABLE ---
add_executable(myfitter myfitter.cc ${sources} ${headers})

# --- 8. LINKING ---
target_link_libraries(myfitter
    PRIVATE
    CHeTLib
    genfit2
    Eigen3::Eigen
    ${ROOT_LIBRARIES}
    ROOT::Geom
    ROOT::Eve
    ${Geant4_LIBRARIES}
    ftxui::screen    # Libreria per gestire il terminale
    ftxui::dom       # Libreria per il layout (box, bordi, testi)
    ftxui::component # Libreria per i widget (bottoni, input, checkbox)
)

# --- 9. SYMLINKS ---
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${PROJECT_SOURCE_DIR}/macros
                ${PROJECT_BINARY_DIR}/macros)
