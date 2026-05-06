cmake_minimum_required(VERSION 3.14)

project(Distiller)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- 1. FIND PACKAGES ---
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree MathCore Physics Hist Matrix Gui Spectrum)
include_directories(${ROOT_INCLUDE_DIRS})

# --- 2. CHET LIBRARY SETUP ---
set(CHeTLibrary_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../CHeTLibrary/install)
include_directories(${CHeTLibrary_DIR}/include)
link_directories(${CHeTLibrary_DIR}/lib)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# --- 3. EXECUTABLE ---
add_executable(Distiller distiller.cc)

# --- 4. COMPILE COMMANDS JSON ---
add_custom_command(TARGET Distiller POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
        ${CMAKE_SOURCE_DIR}/../compile_commands.json
    COMMENT "Copying compile_commands.json to the project root"
)

# --- 5. LINKING ---
target_link_libraries(Distiller
    PRIVATE
    CHeTLib
    ${ROOT_LIBRARIES}
    ROOT::Spectrum
)

# --- 5. EXPORT COMPILE COMMANDS ---
# Copy compile_commands.json to the project root for clangd/LSP support
add_custom_command(TARGET Distiller POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/../compile_commands.json
)
