mesh:
	$(CXX) -c -std=c++11 ../mesh2dh.cpp -Wall
	$(CXX) -c -std=c++11 testmesh.cpp -Wall
	$(CXX) -o testmesh mesh2dh.o testmesh.o -Wall
	
mat:
	$(CXX) -o testmat -std=c++11 testmat.cpp -Wall

run:
	./testmesh
