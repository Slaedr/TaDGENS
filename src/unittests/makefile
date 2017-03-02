mesh:
	$(CXX) -c -std=c++11 ../mesh2dh.cpp
	$(CXX) -c -std=c++11 testmesh.cpp
	$(CXX) -o testmesh mesh2dh.o testmesh.o
	
mat:
	$(CXX) -o testmat -std=c++11 testmat.cpp

run:
	./testmesh
