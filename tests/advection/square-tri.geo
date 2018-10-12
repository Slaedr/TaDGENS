// Simple square test case for linear advection with DGFEM
// Aditya Kashi
// ref to be set from cmd line

refine = 2^ref;
h = 0.25/refine;
xstart = -1.5;
xend = 1.5;
ystart = -1.0;
yend = 1.0;
Point(1) = {xstart, ystart, 0, h};
Point(2) = {xend, ystart, 0, h};
Point(3) = {xend, yend, 0, h};
Point(4) = {xstart, yend, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};
Physical Line(2) = {1, 3};
Physical Line(1) = {2, 4};
Physical Surface(7) = {6};

