//gmsh script
SetFactory("OpenCASCADE");

// parameters
r = DefineNumber[ 0.01, Name "Parameters/r" ]; // radius of the roller
L = DefineNumber[ 0.079, Name "Parameters/L" ]; // height of the elastic gate
S = DefineNumber[ 0.005, Name "Parameters/S" ]; // thickness of the elastic gate
G = DefineNumber[ 0.0025, Name "Parameters/G" ]; // clearance of the elastic gate
Lint = DefineNumber[ 40, Name "Parameters/Lint" ]; // intervals along the length
Sint = DefineNumber[ 8, Name "Parameters/Sint" ]; // intervals along the thickness

Rectangle(1) = {-S, G, 0, S, L, 0}; // X, Y, Z, DX, DY, Rounded radius

Physical Surface("Structure") = {1};
Physical Curve("BeamBottom") = {1};
Physical Curve("BeamRight") = {2};
Physical Curve("BeamTop") = {3};
Physical Curve("BeamLeft") = {4};

Transfinite Curve {1} = Sint;
Transfinite Curve {3} = Sint;
Transfinite Curve {2} = Lint;
Transfinite Curve {4} = Lint;

Transfinite Surface {1};
Recombine Surface{1};
