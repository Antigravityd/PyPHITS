from angel import *
from parsec import *


print(inps.parse("INPS: { file.eps} X( 10 ) Y( 11 ) IX( 12.3 ) IY( 55.3 ) S(1.2) A(30 )"))

print(inps.parse("INPS: { file.eps} X( 10 ) Y( 11 ) IX( 12.3 ) IY( 55.3 ) SX(1.) SY( 2.5 ) A(30 )"))

print(set_stmt.parse(" SET: c1[ 10*2 ]"))

print(x.parse("X: Title X"))

print(y.parse("Y: Title Y :)"))

print(normal_h.parse("""'Graph'
X: Xtitle
Y: Ytitle
H: X  Y1*10+12 NY2(Error) D1
   1  2        3   4
   5  6        7   8
   9  10       11  12
"""))

print(normal_h.parse("""'Graph'
X: Xtitle
Y: Ytitle
H: X  Y1(Why One),C[H(0.5)S(0.5)B(0.5)]12[C[darkred]]XS[12]HHH NY2 D1**3
   1  2        3   4
   5  6        7   8
   9  10       11  12
"""))

print(contour.parse("""'Graph'
X: Xtitle
Y: Ytitle
H2: Y = 3.0 TO -3.0 BY -1.0 ; X = -3.0 TO 3.0 BY 1.0 ;
0 0 0 0 0 0 0
0 8 0 2 4 0 0
0 0 1 5 6 2 0
0 2 4 14 15 9 0
0 2 5 23 32 12 0
0 0 0 0 0 0 0"""))

print(cluster.parse("""'Graph'
X: Xtitle
Y: Ytitle
HD: Y = 3.0 TO -3.0 BY -1.0 ; X = -3.0 TO 3.0 BY 1.0 ;
0 0 0 0 0 0 0
0 8 0 2 4 0 0
0 0 1 5 6 2 0
0 2 4 14 15 9 0
0 2 5 23 32 12 0
0 0 0 0 0 0 0"""))

print(color_cluster.parse("""'Graph'
X: Xtitle
Y: Ytitle
HC: Y = 3.0 TO -3.0 BY -1.0 ; X = -3.0 TO 3.0 BY 1.0 ;
0 0 0 0 0 0 0
0 8 0 2 4 0 0
0 0 1 5 6 2 0
0 2 4 14 15 9 0
0 2 5 23 32 12 0
0 0 0 0 0 0 0"""))

print(comment.parse("W:Comment/ X(10) Y( 12) IX(2 ) IY(2) S(1.1) C(cyan) a(12.3)"))


print(multiline_comment.parse("""WT: X(10) Y( 12) IX(2 ) IY(2) S(1.1) C(cyan) A(12.3) B(12)
this
is a
multiline
comment.
E:
"""))

print(tabular_comment.parse("""WTAB: TAB{ 6 } X(10) Y( 12) IX(2 ) IY(2) S(1.1) C(cyan) A(12.3) CB(cyan)
this
is a
multiline
comment.
E:
"""))

print(comment_arrow.parse("AW: comment arrow / X(10) Y( 12) AX(12.3) AY(23) IL IX(2 ) IY(2) S(1.1) C(cyan) TTT"))

print(arrow.parse("A: X(12.3) Y(22.434) AX(2343) AY(1.200) C(blue) A(121) TT"))

print(hollow_arrow.parse("AB: X(20.3) Y(12.5) AX(2121) AY(1.200) C(red) A(12) N"))

print(polygon.parse("POLG: X(20.3) Y(12.5) SX(2121) SY(   1.200) CL(red) CB(green) A(12) PL(6)"))

print(box.parse("BOX: X(20.3) Y(12.5) S(123) CL(red) CB(green) CS(blue) A(12) BOX(bruh)"))

print(ribbon.parse("RIBN: X(20.3) Y(12.5) S(123) CL(red) CB(green) CS(blue) A(12) ZZZ"))

print(star.parse("STAR: X(20.3) Y(12.5) S(123) CL(red) CB(green) CS(blue) A(12) ZZZ PL(5) V(12)"))

print(parameter.parse("P: FRAM(1)"))

print(page.parse("""'Figure 2'
X: E(MeV)
Y: d\sigma/d\OmegadE (nb/MeV/sr)
W:Comment=^{12}C + ^{40}Ca/ X(1) Y(65)"""))
# A: X(4.3) Y(65) AX(5.5) AY(53)
# AW:Arrow/ X(4.5) Y(23) AX(4) AY(10)
# H: X Y(Histogram),DH0
# 0.0 0
# 0.5 10
# 1.5 20
# 2.5 35
# 3.5 30
# 4.5 45
# 5.5 70
# 6.5 0
# Test comment
# H: X Y(Spline),SL5
# 0.5 40
# 4 10
# 6 20
# 7 60
# H: X DX Y(Error bar),N3 D
# 2 1 50 10"""))
