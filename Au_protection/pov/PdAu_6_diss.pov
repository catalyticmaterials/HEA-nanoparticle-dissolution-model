#version 3.6;
#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 2.2 max_trace_level 6}
background {color Grey transmit 1.0}
camera {orthographic
  right -10.98*x up 10.98*y
  direction 1.00*z
  location <0,0,50.00> look_at <0,0,0>}


light_source {<  -30.00,  30.00,   40.00> color Gray40 shadowless}
light_source {<  30.00,  30.00,   40.00> color Gray40 shadowless} 
light_source {<  30.0,  -30.00,   40.00> color Gray40 shadowless}
light_source {<  -30.0,  -30.00,   40.00> color Gray40 shadowless} 
light_source {<  0.0,  0.00,   40.00> color Gray25 shadowless}
// no fog
#declare simple = finish {phong 0.7 ambient 0.4 diffuse 0.55}
#declare pale = finish {ambient 0.9 diffuse 0.30 roughness 0.001 specular 0.2 }
#declare intermediate = finish {ambient 0.4 diffuse 0.6 specular 0.1 roughness 0.04}
#declare vmd = finish {ambient 0.2 diffuse 0.80 phong 0.25 phong_size 10.0 specular 0.2 roughness 0.1}
#declare jmol = finish {ambient 0.4 diffuse 0.6 specular 1 roughness 0.001 metallic}
#declare ase2 = finish {ambient 0.2 brilliance 3 diffuse 0.6 metallic specular 0.7 roughness 0.04 reflection 0.15}
#declare ase3 = finish {ambient 0.4 brilliance 2 diffuse 0.6 metallic specular 1.0 roughness 0.001 reflection 0.0}
#declare glass = finish {ambient 0.4 diffuse 0.35 specular 1.0 roughness 0.001}
#declare glass2 = finish {ambient 0.3 diffuse 0.3 specular 1.0 reflection 0.25 roughness 0.001}
#declare Rcell = 0.100;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
     torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
     translate LOC}
#end

// no cell vertices
atom(< -1.93,  -1.82,  -4.57>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #0
atom(< -0.40,   1.53,  -1.49>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #1
atom(< -1.90,   2.02,  -3.77>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #2
atom(< -4.26,  -2.10,  -3.10>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #3
atom(< -1.50,   0.49,  -6.04>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #4
atom(< -3.83,   0.20,  -4.57>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #5
atom(< -2.33,  -0.29,  -2.29>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #6
atom(<  0.81,  -3.07,  -8.31>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #7
atom(< -0.44,  -2.31,  -2.28>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #8
atom(<  2.74,  -1.25,  -7.50>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #9
atom(< -0.03,  -3.84,  -4.56>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #10
atom(<  0.40,  -1.53,  -6.03>, 1.36, rgbt <0.25, 0.41, 0.88, 0.00>, 0.0, ase3) // #11
atom(<  1.90,  -2.02,  -3.75>, 1.36, rgbt <0.25, 0.41, 0.88, 0.00>, 0.0, ase3) // #12
atom(<  2.33,   0.29,  -5.23>, 1.36, rgbt <0.25, 0.41, 0.88, 0.00>, 0.0, ase3) // #13
atom(<  0.00,   0.00,  -3.76>, 1.36, rgbt <0.25, 0.41, 0.88, 0.00>, 0.0, ase3) // #14
atom(< -0.84,  -0.77,  -0.01>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #15
atom(<  1.93,   1.82,  -2.95>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #16
atom(<  0.44,   2.31,  -5.24>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #17
atom(<  4.23,  -1.73,  -5.22>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #18
atom(<  2.30,  -3.55,  -6.02>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #19
atom(<  1.06,  -2.79,   0.00>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #20
atom(<  3.83,  -0.20,  -2.94>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #21
atom(<  1.50,  -0.49,  -1.48>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #22
atom(<  0.84,   0.77,  -7.51>, 1.36, rgbt <1.00, 0.84, 0.00, 0.00>, 0.0, ase3) // #23

// no constraints
