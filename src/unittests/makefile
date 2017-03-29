
CXXFLAGS := -std=c++14 -Wall -ggdb -DDEBUG -I${EIGEN_DIR}

.PHONY: run

amesh2dh.o: ../amesh2dh.cpp
	$(CXX) -c ${CXXFLAGS} ../amesh2dh.cpp

aquadrature.o: ../aquadrature.cpp
	$(CXX) -c ${CXXFLAGS} ../aquadrature.cpp

aelements.o: ../aelements.cpp
	$(CXX) -c ${CXXFLAGS} ../aelements.cpp
	
mat: testmat.cpp
	$(CXX) -o mat -std=c++14 testmat.cpp -Wall

mesh: amesh2dh.o testmesh.cpp
	$(CXX) -c ${CXXFLAGS} testmesh.cpp
	$(CXX) -o mesh amesh2dh.o testmesh.o

element: aelements.o aquadrature.o amesh2dh.o testelement.cpp
	$(CXX) -c ${CXXFLAGS} testelement.cpp
	$(CXX) -o element aquadrature.o aelements.o amesh2dh.o testelement.o

run:
	./testmesh
