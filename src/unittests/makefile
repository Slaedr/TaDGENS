mesh:
	$(CC) -c -std=c++11 ../mesh2dh.cpp
	$(CC) -c -std=c++11 testmesh.cpp
	$(CC) -o testmesh mesh2dh.o testmesh.o
	
mat:
	$(CC) -o testmat -std=c++11 testmat.cpp

run:
	./testmesh