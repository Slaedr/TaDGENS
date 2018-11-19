// Geometry for investigating curved outflow boundary for linear advection with DGFEM
// Aditya Kashi
// ref to be set from cmd line

refine = 2^ref;
h = 0.5/refine;
xstart = -1.5;
xend = 1.5;
ystart = -1.0;
yend = 1.0;

Point(1) = {xstart, ystart, 0, h};
Point(2) = {xend, ystart, 0, h};
Point(3) = {xend, yend, 0, h};
Point(4) = {xstart, yend, 0, h};
Point(5) = {xend, (ystart+yend)/2.0, 0, h};
Point(6) = {xend+(xend-xstart)*0.4, (ystart+yend)/2.0, 0, h};

Line(1) = {1, 2};
Ellipse(2) = {2, 5, 2, 6};
Ellipse(3) = {6, 5, 2, 3};
Line(4) = {3, 4};
Line(5) = {4, 1};

Line Loop(6) = {1, 2, 3, 4, 5};
Plane Surface(6) = {6};
Physical Line(2) = {1, 4};
Physical Line(1) = {2, 3, 5};
Physical Surface(7) = {6};

Color Black { Surface {6}; }
