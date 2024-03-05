// gmsh script
//+
SetFactory("OpenCASCADE");
//+
Circle(1) = {0.2, 0.2, 0, 0.05, -Pi, Pi};
Point(2) = {0.2, 0.21, 0};
Point(3) = {0.6, 0.21, 0};
Point(4) = {0.6, 0.19, 0};
Point(5) = {0.2, 0.19, 0};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
//+
BooleanFragments{ Curve{1}; Delete; }{ Line{2}; Line{4}; Delete; }
Recursive Delete {
  Curve{4}; Curve{6}; Curve{7}; Curve{10}; 
}
//+
Curve Loop(1) = {5, 8, 3, 9};
Surface(1) = {1};
//+
Physical Surface("FLAG") = {1};
Physical Curve("FixedDisplacement") = {5};
Physical Curve("FlagTop") = {8};
Physical Curve("FlagBottom") = {9};
Physical Curve("FlagEnd") = {3};
//+
Transfinite Curve {3, 5} = 15;
Transfinite Curve {8, 9} = 175;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};

