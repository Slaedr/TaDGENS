// Simple square test case for linear advection with DGFEM
// Aditya Kashi

refine = 2^ref;
h = 0.25/refine;
npoin = 10*refine;
xstart = -1.5;
xend = 1.5;
ystart = -1.0;
yend = 1.0;
Point(1) = {xstart, ystart, 0, h};
Point(2) = {xend, ystart, 0, h};
Point(3) = {xend, yend, 0, h};
Point(4) = {xstart, yend, 0, h};
Line(1) = {1, 2}; Transfinite Line(1) = npoin;
Line(2) = {2, 3}; Transfinite Line(2) = npoin;
Line(3) = {3, 4}; Transfinite Line(3) = npoin;
Line(4) = {4, 1}; Transfinite Line(4) = npoin;
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6}; Transfinite Surface{6};
Physical Line(2) = {1, 3};
Physical Line(1) = {2, 4};
Physical Surface(7) = {6};

// Comment out for triangular mesh
Recombine Surface {6};
