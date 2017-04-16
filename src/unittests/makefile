
CXXFLAGS := -std=c++14 -Wall -ggdb -DDEBUG -I${EIGEN_DIR}

.PHONY: run
.PHONY: clean

amesh2dh.o: ../amesh2dh.cpp
	${CXX} -c ${CXXFLAGS} ../amesh2dh.cpp

aquadrature.o: ../aquadrature.cpp
	${CXX} -c ${CXXFLAGS} ../aquadrature.cpp

aelements.o: ../aelements.cpp
	${CXX} -c ${CXXFLAGS} ../aelements.cpp
	
mat: testmat.cpp
	${CXX} ${CXXFLAGS} -o mat testmat.cpp

mesh: amesh2dh.o testmesh.cpp
	${CXX} -c ${CXXFLAGS} testmesh.cpp
	${CXX} -o mesh amesh2dh.o testmesh.o

elementtri: aelements.o aquadrature.o amesh2dh.o testelementtri.cpp
	${CXX} -c ${CXXFLAGS} testelementtri.cpp
	${CXX} -o elementtri aquadrature.o aelements.o amesh2dh.o testelementtri.o

run:
	./mesh
	./elementtri

clean:
	rm *.o
	rm mesh
	rm elementtri
