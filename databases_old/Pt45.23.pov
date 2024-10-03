#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {orthographic
  right -17.31*x up 17.31*y
  direction 1.00*z
  location <0,0,50.00> look_at <0,0,0>}


light_source {<  2.00,   3.00,  40.00> color White
  area_light <0.70, 0, 0>, <0, 0.70, 0>, 3, 3
  adaptive 1 jitter}
// no fog
#declare simple = finish {phong 0.7}
#declare pale = finish {ambient 0.5 diffuse 0.85 roughness 0.001 specular 0.200 }
#declare intermediate = finish {ambient 0.3 diffuse 0.6 specular 0.1 roughness 0.04}
#declare vmd = finish {ambient 0.0 diffuse 0.65 phong 0.1 phong_size 40.0 specular 0.5 }
#declare jmol = finish {ambient 0.2 diffuse 0.6 specular 1 roughness 0.001 metallic}
#declare ase2 = finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic specular 0.7 roughness 0.04 reflection 0.15}
#declare ase3 = finish {ambient 0.15 brilliance 2 diffuse 0.6 metallic specular 1.0 roughness 0.001 reflection 0.0}
#declare glass = finish {ambient 0.05 diffuse 0.3 specular 1.0 roughness 0.001}
#declare glass2 = finish {ambient 0.01 diffuse 0.3 specular 1.0 reflection 0.25 roughness 0.001}
#declare Rcell = 0.020;
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
atom(< -3.89,  -4.95,  -6.03>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #0
atom(< -1.06,  -4.95,  -6.03>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #1
atom(<  1.76,  -4.95,  -6.03>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #2
atom(< -2.48,  -3.72,  -8.14>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #3
atom(<  0.35,  -3.72,  -8.14>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #4
atom(<  3.17,  -3.72,  -8.14>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #5
atom(< -1.06,  -2.50, -10.26>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #6
atom(<  1.76,  -2.50, -10.26>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #7
atom(<  4.59,  -2.50, -10.26>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #8
atom(< -3.89,  -3.76,  -3.46>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #9
atom(< -1.06,  -3.76,  -3.46>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #10
atom(<  1.76,  -3.76,  -3.46>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #11
atom(< -2.48,  -2.54,  -5.58>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #12
atom(<  0.35,  -2.54,  -5.58>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #13
atom(<  3.17,  -2.54,  -5.58>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #14
atom(< -1.06,  -1.32,  -7.70>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #15
atom(<  1.76,  -1.32,  -7.70>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #16
atom(<  4.59,  -1.32,  -7.70>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #17
atom(< -2.48,  -1.36,  -3.01>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #18
atom(<  0.35,  -1.36,  -3.01>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #19
atom(<  3.17,  -1.36,  -3.01>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #20
atom(< -1.06,  -0.14,  -5.13>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #21
atom(<  1.76,  -0.14,  -5.13>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #22
atom(<  4.59,  -0.14,  -5.13>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #23
atom(<  0.35,   1.09,  -7.25>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #24
atom(<  3.17,   1.09,  -7.25>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #25
atom(<  6.00,   1.09,  -7.25>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #26
atom(< -3.89,   1.05,  -2.57>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #27
atom(< -1.06,   1.05,  -2.57>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #28
atom(<  1.76,   1.05,  -2.57>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #29
atom(< -2.48,   2.27,  -4.68>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #30
atom(<  0.35,   2.27,  -4.68>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #31
atom(<  3.17,   2.27,  -4.68>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #32
atom(< -1.06,   3.49,  -6.80>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #33
atom(<  1.76,   3.49,  -6.80>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #34
atom(<  4.59,   3.49,  -6.80>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #35
atom(< -3.89,   2.23,   0.00>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #36
atom(< -1.06,   2.23,   0.00>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #37
atom(<  1.76,   2.23,   0.00>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #38
atom(< -2.48,   3.45,  -2.12>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #39
atom(<  3.17,   3.45,  -2.12>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #40
atom(< -1.06,   4.67,  -4.24>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #41
atom(<  1.76,   4.67,  -4.24>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #42
atom(<  4.59,   4.67,  -4.24>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, jmol) // #43

// no constraints
