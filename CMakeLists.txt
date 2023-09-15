cmake_minimum_required(VERSION 3.12)

project (LBM_cluster)

set(PROJECT_NAME LBM_cluster)

set(INCLUDES
	includes/Solvers/BasicSolver2D.hpp
	includes/Solvers/SolverMRTDimensional2D.hpp
	includes/Solvers/Solver3Phases2D.hpp
	includes/Solvers/SolverPVTsim2D.hpp
	includes/Solvers/SolverPVTsimDynamic2D.hpp
)

set(SOURCE_EXE
	src/main.cpp
)

set(SOURCES_SOLVERS
	src/Solvers/BasicSolver2D.cpp
	src/Solvers/SolverMRTDimensional2D.cpp
	src/Solvers/Solver3Phases2D.cpp
	src/Solvers/SolverPVTsim2D.cpp
	src/Solvers/SolverPVTsimDynamic2D.cpp
)

set(ALL_SOURCES
	${INCLUDES}
	${SOURCES_SOLVERS}
)

add_library(LBM_core STATIC
${ALL_SOURCES}
)



include_directories(${PROJECT_NAME} PUBLIC includes)
include_directories(${PROJECT_NAME} PRIVATE src)
include_directories(${PROJECT_NAME} PRIVATE VTK)
include_directories(${PROJECT_NAME} PRIVATE Masks)
set (CMAKE_CXX_STANDARD 17)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES 
	${ALL_SOURCES}
)

add_executable(main ${SOURCE_EXE})