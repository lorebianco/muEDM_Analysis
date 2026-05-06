cmake_minimum_required(VERSION 3.16)

project(ToyMC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- 1. FIND PACKAGES ---
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist MathCore Tree ROOTDataFrame)
include_directories(${ROOT_INCLUDE_DIRS})

# --- 1.1 FETCHCONTENT: FTXUI ---
include(FetchContent)
FetchContent_Declare(ftxui
  GIT_REPOSITORY https://github.com/ArthurSonzogni/FTXUI
  GIT_TAG v5.0.0
)
FetchContent_MakeAvailable(ftxui)

# --- 2. CHET LIBRARY SETUP ---
set(CHeTLibrary_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../CHeTLibrary/install)
include_directories(${CHeTLibrary_DIR}/include)
link_directories(${CHeTLibrary_DIR}/lib)

# --- 3. PROJECT INCLUDES ---
include_directories(${PROJECT_SOURCE_DIR}/include)

# --- 4. SOURCES ---
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# --- 5. EXECUTABLE ---
# Assumes the main entry point will be toymc.cc in the ToyMC directory
add_executable(toymc toymc.cc ${sources} ${headers})

# --- 6. LINKING ---
target_link_libraries(toymc
    PRIVATE
    CHeTLib
    ${ROOT_LIBRARIES}
    ftxui::screen
    ftxui::dom
    ftxui::component
)
