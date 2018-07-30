cmake_minimum_required(VERSION 2.6)
project(generalhbmcontinuation)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_BUILD_TYPE Release)

find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)
find_package(PythonLibs REQUIRED)


set(TRILINOS_ROOT /home/jiri/Libraries/Trilinos_pytrilinos/)
set(TRILINOS_BUILD_PACKAGES /home/jiri/Libraries/Trilinos_pytrilinos/build/packages/)
set(PETSC_BASE /home/jiri/Libraries/PETSc/petsc-3.9.2/)

include_directories(${TRILINOS_ROOT}packages/amesos/src)
include_directories(${TRILINOS_ROOT}packages/epetra/src)
include_directories(${TRILINOS_ROOT}build/packages/epetra/src)
include_directories(${TRILINOS_ROOT}packages/nox/src)
include_directories(${TRILINOS_ROOT}packages/nox/src-epetra)
include_directories(${TRILINOS_ROOT}packages/nox/src-lapack)
include_directories(${TRILINOS_ROOT}build/packages/nox/src)
include_directories(${TRILINOS_ROOT}build/packages/nox/src-lapack)
include_directories(${TRILINOS_ROOT}packages/teuchos/core/src)
include_directories(${TRILINOS_ROOT}packages/teuchos/comm/src)
include_directories(${TRILINOS_ROOT}packages/teuchos/numerics/src)
include_directories(${TRILINOS_ROOT}packages/teuchos/parameterlist/src)
include_directories(${TRILINOS_ROOT}build/packages/teuchos/core/src)
include_directories(${TRILINOS_ROOT}packages/ml/src/Main)
include_directories(${TRILINOS_ROOT}packages/ml/src/Utils)
include_directories(${TRILINOS_ROOT}packages/ml/src/Include)
include_directories(${TRILINOS_ROOT}packages/ml/src/Smoother)
include_directories(${TRILINOS_ROOT}packages/ml/src/Comm)
include_directories(${TRILINOS_ROOT}packages/ml/src/Operator)
include_directories(${TRILINOS_ROOT}packages/ml/src/Krylov)
include_directories(${TRILINOS_ROOT}packages/ml/src/Coarsen)
include_directories(${TRILINOS_ROOT}packages/ml/src/FEGrid)
include_directories(${TRILINOS_ROOT}build/packages/ml/src)
include_directories(${TRILINOS_ROOT}packages/aztecoo/src)
include_directories(${TRILINOS_ROOT}build/packages/aztecoo/src)
include_directories(${TRILINOS_ROOT}packages/epetraext/src)
include_directories(${TRILINOS_ROOT}packages/epetraext/src/transform)
include_directories(${TRILINOS_ROOT}build/packages/epetraext/src)
include_directories(${TRILINOS_ROOT}packages/anasazi/src)
include_directories(${TRILINOS_ROOT}packages/anasazi/epetra/src)
include_directories(${TRILINOS_ROOT}build/packages/anasazi/epetra/src)
include_directories(${TRILINOS_ROOT}build/packages/anasazi/src)
include_directories(${TRILINOS_ROOT}packages/nox/src-loca/src)
include_directories(${TRILINOS_ROOT}packages/nox/src-loca/src-lapack)

include_directories(${MPI_INCLUDE_PATH})
include_directories(${PYTHON_INCLUDE_DIRS})

include_directories(${PETSC_BASE}include)
include_directories(${PETSC_BASE}arch-linux2-c-debug/include)
include_directories("/home/jiri/Libraries/matplotlib-cpp-master")


link_directories(${TRILINOS_BUILD_PACKAGES}nox/src)
link_directories(${TRILINOS_BUILD_PACKAGES}nox/src-epetra)
link_directories(${TRILINOS_BUILD_PACKAGES}nox/src-lapack)
link_directories(${TRILINOS_BUILD_PACKAGES}epetra/src)
link_directories(${TRILINOS_BUILD_PACKAGES}epetraext/src)
link_directories(${TRILINOS_BUILD_PACKAGES}teuchos/core/src)
link_directories(${TRILINOS_BUILD_PACKAGES}teuchos/comm/src)
link_directories(${TRILINOS_BUILD_PACKAGES}teuchos/numerics/src)
link_directories(${TRILINOS_BUILD_PACKAGES}teuchos/parameterlist/src)
link_directories(${TRILINOS_BUILD_PACKAGES}teuchos/remainder/src)
link_directories(${TRILINOS_BUILD_PACKAGES}aztecoo/src)
link_directories(${TRILINOS_BUILD_PACKAGES}ml/src)
link_directories(${TRILINOS_BUILD_PACKAGES}thyra/core/src)
link_directories(${TRILINOS_BUILD_PACKAGES}thyra/adapters/epetra/src)
link_directories(${TRILINOS_BUILD_PACKAGES}thyra/adapters/epetraext/src)
link_directories(${TRILINOS_BUILD_PACKAGES}ifpack/src)
link_directories(${TRILINOS_BUILD_PACKAGES}kokkos/core/src)
link_directories(${TRILINOS_BUILD_PACKAGES}zoltan/src)
link_directories(${TRILINOS_BUILD_PACKAGES}amesos/src)
link_directories(${TRILINOS_BUILD_PACKAGES}common/auxiliarySoftware/SuiteSparse/src)
link_directories(${TRILINOS_BUILD_PACKAGES}anasazi/src)
link_directories(${TRILINOS_BUILD_PACKAGES}nox/src-loca/src)
link_directories(${TRILINOS_BUILD_PACKAGES}nox/src-loca/src-lapack)

add_executable(generalhbmcontinuation main.cpp Functions.cpp ProblemInterface.cpp Nonlinearities/NonlinearBase.cpp Nonlinearities/CubicSpring.cpp Nonlinearities/NonlinearBaseFD.cpp ContinuationWrapper.cpp Nonlinearities/FricElem3D.cpp Nonlinearities/PenaltyWall.cpp TestingFunctions.cpp Time.cpp Config.cpp AftSimple.cpp)

set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HBMCont.exe")

target_link_libraries(generalhbmcontinuation ${MPI_LIBRARIES} 
	ml zoltan loca localapack nox noxepetra noxlapack epetra ml aztecoo
	epetraext teuchoscore 
	teuchoscomm teuchosnumerics teuchosparameterlist 
	thyracore thyraepetra thyraepetraext amesos trilinosss
	kokkoscore dl fftw3 fftw3_omp
	ifpack teuchosremainder ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${PYTHON_LIBRARIES} -lmesh_devel)


set_target_properties(generalhbmcontinuation PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")

set_target_properties(generalhbmcontinuation PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

