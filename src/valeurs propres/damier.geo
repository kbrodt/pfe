h = 0.2;
a = 1;
Point(1) = {0, 0, 0, 0.05*h};
Point(2) = {a, 0, 0, h};
Point(3) = {a, a, 0, h};
Point(4) = {0, a, 0, h};
Point(5) = {-a, a, 0, h};
Point(6) = {-a, 0, 0, h};
Point(7) = {-a, -a, 0, h};
Point(8) = {0, -a, 0, h};
Point(9) = {a, -a, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {6, 1};
Line(6) = {4, 5};
Line(7) = {5, 6};
Line(8) = {7, 8};
Line(9) = {8, 1};
Line(10) = {6, 7};
Line(11) = {8, 9};
Line(12) = {9, 2};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, -4, 6, 7};
Line Loop(3) = {8, 9, -5, 10};
Line Loop(4) = {11, 12, -1, -9};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

Physical Point(1) = {2,3,4,5,6,7,8,9};

Physical Line(1) = {2, 3, 6, 7, 10, 8, 11, 12};

Physical Surface(1) = {1, 3};
Physical Surface(2) = {2, 4};
