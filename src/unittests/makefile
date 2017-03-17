EIGEN_PATH :=  ../../eigen-3.3.3
CXXFLAGS := -std=c++14 -Wall -ggdb -I${EIGEN_PATH} -DDEBUG

meshlib:
	$(CXX) -c ${CXXFLAGS} ../amesh2dh.cpp
	
mat:
	$(CXX) -o mat -std=c++14 testmat.cpp -Wall

quadraturelib:
	$(CXX) -c ${CXXFLAGS} ../aquadrature.cpp

elemlib:
	$(CXX) -c ${CXXFLAGS} ../aelements.cpp

mesh: meshlib
	$(CXX) -c ${CXXFLAGS} testmesh.cpp
	$(CXX) -o mesh amesh2dh.o testmesh.o

element: elemlib quadraturelib meshlib
	$(CXX) -c ${CXXFLAGS} testelement.cpp
	$(CXX) -o element aquadrature.o aelements.o amesh2dh.o testelement.o

run:
	./testmesh
