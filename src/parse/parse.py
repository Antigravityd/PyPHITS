from lark import Lark

test = r"""[Title]
2023-02-22 09:26:52.596846
[Parameters]
negs = 1
e-mode = 2
[Source]
s-type = 1
nproj = 241Am
ispfs = 2
z0 = -0.25
z1 = 0.25
r0 = 0.3
e0 = 2.2
[Material]
MAT[1]
C 0.65 H 0.2333333333333333 O 0.11666666666666665

MAT[2]
C 0.64 H 0.24 O 0.12

MAT[3]
C 0.56 H 0.29333333333333333 O 0.14666666666666667

MAT[4]
C 0.57 H 0.2866666666666666 O 0.1433333333333333

mat[5]
C 0.62 H 0.2533333333333333 O 0.12666666666666665

MAT[6]
C 0.63 H 0.24666666666666665 O 0.12333333333333332

MAT[7]
C 0.61 H 0.26 O 0.13

MAT[8]
C 0.55 H 0.3 O 0.15

MAT[9]
C 0.59 H 0.2733333333333333 O 0.13666666666666666

MAT[10]
C 0.58 H 0.28 O 0.14

MAT[11]
C 0.6 H 0.26666666666666666 O 0.13333333333333333

[Surface]
1 SO 2
2 SO 5
3 SO 10
4 SO 6
5 SO 4
6 SO 7
7 SO 11
8 SO 9
9 SO 3
10 SO 8
11 SO 1
[Cell]
1 7 -1  6 4


2 2 -1  3 8


3 8 -1  11


4 9 -1  2 5


5 -1  7


6 11 -1  4 2


7 5 -1  10 6


8 1 -1  7 3


9 4 -1  9 1


10 10 -1  5 9


11 6 -1  8 10


12 3 -1  1 11


[T-Cross]
mesh = reg
unit = 1
axis = reg
file = cross.dmp
output_type = flux
reg = 1
dump = -2
1 8
mesh = reg
unit = 1
axis = reg
file = cross.dmp
output_type = flux
reg = 1
dump = -2
1 8
r-out r-in area
8
5
area = 4.1887902047863905
$ This could be e.g. an [Elastic Option] section text.
"""
parser = Lark(r"""
// PREPROCESSING EXPECTATIONS: comments removed, only spaces and newlines as whitespace, divided lines expanded,
// continued lines continued, bad section lines and off sections removed, terminators, skipped lines,
// files inserted, user-defined variables evaluated and substituted, mathematical expressions evaluated and substituted,
// and blank lines removed.
// TODO: revise this to ensure only the truly significantly context-dependent parts are removed.
// See the 'reading control' portion of the PHITS manual.
start: section*

section: title
       | parameters
       | sources
       | materials
       | surfaces
       | cells
       | transforms
       // | temperature
       // | mattimechange
       // | magneticfield
       // | electromagneticfield
       // | deltaray
       // | trackstructure
       // | supermirror
       // | elasticoption
       // | datamax
       // | fragdata
       // | importance
       // | weightwindow
       // | wwbias
       // | forcedcollisions
       // | repeatedcollisions
       // | volume
       // | multiplier
       // | matnamecolor
       // | regname
       // | counter
       // | timer
       // | ttrack
       // | tcross
       // | tpoint
       // | tdeposit
       // | tdeposittwo
       // | theat
       // | tyield
       // | tproduct
       // | tdpa
       // | tlet
       // | tsed
       // | ttime
       // | tinteract
       // | tdchain
       // | twwg
       // | twwbg
       // | tvolume
       // | tuserdefined
       // | tgshow
       // | trshow
       // | tthreedshow


title: "[" "t" "i" "t" "l" "e" "]" /[^[]+/

parameters: "[" "p" "a" "r" "a" "m" "e" "t" "e" "r" "s" "]" assignment*




sources: "[" "s" "o" "u" "r" "c" "e" "]" [source
                                          | (assignment* ("<source>" "=" computation source)+ assignment*)]




source: "s-type" "=" POSINT (assignment | etype | ttype | atype)*



materials: "[" "m" "a" "t" "e" "r" "i" "a" "l" "]" columntitle? (material | thermlib | chemical)+

material: ("MAT[" POSINT "]" |"M" POSINT) ((element computation) | assignmentinl)+

thermlib: "MT" POSINT /[\.0-9a-z]/i
element: (/[0-9]*[a-z]+/ | /[a-z]+-[0-9]+/ | /[0-9]+/) /\.[0-9][0-9][a-z]/i
chemical: "chem" "="  ((element computation) | assignmentinl)+



_SURFNAME: / *\[ *s *u *r *f *a *c *e *\] *\n/i
surfaces: "[" "s" "u" "r" "f" "a" "c" "e" "]" (((computation | surfsymbol)+ )+

surfsymbol: (/[a-z]/ | "/")+


// lattice-universe is sus, as always
_CELLNAME: / *\[ *c *e *l *l *\] *\n/i
cells: _CELLNAME (" "* ((computation " "+) ~ 3 celldef | computation " "+ likebut) " "+ closeassignmentinl " "* "\n")+

celldef: INT | celldef " "+ celldef | celldef " "* ":" " "* celldef | "#" celldef | "(" " "* celldef " "* ")"

likebut: "like"i " "+ POSINT " "+ "but"i " "+


_TFNAME: / *\[ *t *r *a *n *s *f *o *r *m *\] *\n/i
transforms: _TFNAME " "* "*"? "tr"i POSINT (" " | "\n")+ (computation (" " | "\n")+) ~ 13



assignment: " "* iden / *= */ vals / *\n/
assignmentinl: iden / *= */ vals
closeassignmentinl: iden "=" vals
iden: IDENTIFIER | (IDENTIFIER "(" POSINT ")")
vals: computation | FILENAME

regiongrp: POSINT | /( */ (POSINT " "*) / *)/ | /\{ */ POSINT
latuniv: // TODO; lattice universe specification like ( 6 < 10[1 0 0] < u=3 )

etype: e1 | e2 | e3 | e4 | e5 | e6 | e7 | e20 | e25 | e28

e1: " "* "e-type"i " "* "=" " "* /1|8|11|18|21|22|31|32/ " "* "\n" assignment grid

e4: " "* "e-type"i " "* "=" " "* /4|9|14|19|23|24|33|34/  " "* "\n" assignment grid assignment grid

e2: " "* "e-type"i " "* "=" " "* /2|12/  " "* "\n" assignment ~ 4

e3: " "* "e-type"i " "* "=" " "* /3/  " "* "\n" assignment ~ 5

e5: " "* "e-type"i " "* "=" " "* /5|15/  " "* "\n" assignment ~ 4


e6: " "* "e-type"i " "* "=" " "* /6|16/  " "* "\n" assignment ~ 5 grid

e7: " "* "e-type"i " "* "=" " "* /7/  " "* "\n" assignment ~ 6 grid

e20: " "* "e-type"i " "* "=" " "* /20/  " "* "\n" assignment

e25: " "* "e-type"i " "* "=" " "* /25|26/  " "* "\n" assignment ~ 14

e28: " "* "e-type"i " "* "=" " "* /28|29/  " "* "\n" assignment ~ 7



atype: p1 | p4 | p5 | p6

p1: " "* "a-type"i " "* "=" " "* /1|11/  " "* "\n" assignment grid


p4: " "* "a-type"i " "* "=" " "* /4|14/  " "* "\n" assignment grid assignment grid

p5: " "* "a-type"i " "* "=" " "* /5|15/  " "* "\n" assignment ~ 4

p6: " "* "a-type"i " "* "=" " "* /6|16/  " "* "\n" assignment ~ 3 grid


ttype: t0 | t3 | t4 | t5 | t6  | t100

t0: " "* "t-type"i " "* "=" " "* /0|1|2/  " "* "\n" assignment ~ 5


t3: " "* "t-type"i " "* "=" " "* /3/  " "* "\n" assignment grid

t4: " "* "t-type"i " "* "=" " "* /0|1|2/  " "* "\n" assignment grid assignment grid

t5: " "* "t-type"i " "* "=" " "* /5/  " "* "\n" assignment ~ 4

t6: " "* "t-type"i " "* "=" " "* /6/  " "* "\n" assignment ~ 5 grid

t100: " "* "t-type"i " "* "=" " "* /100/  " "* "\n" assignment ~ 2




grid: ((computation " "+)+ "\n")+
columntitle: (IDENTIFIER " "+)+

IDENTIFIER: /[a-z0-9-]+/i
POSINT: /[0-9]+/
INT: /-?[0-9]+/

function{name}: /name/i "(" computation ")"
binop{sym}: computation sym computation

computation: " "* (NUMBER | /pi/i | /x/i | USRCONSTANT
  	       	  | binop{/\+/} | binop{/-/} | binop{/\*/} | binop{/\//} | binop{/\*\*/}
	      	  | function{/float/i} | function{/int/i} | function{/abs/i} | function{/exp/i} | function{/log/i} | function{/log10/i}
	      	  | function{/max/i} | function{/min/i} | function{/mod/i} | function{/nint/i} | function{/sign/i} | function{/sqrt/i}
	      	  | function{/acos/i} | function{/asin/i} | function{/atan/i} | function{/atan2/i} | function{/cos/i}
	      	  | function{/cosh/i} | function{/sin/i} | function{/sinh/i} | function{/tan/i} | function{/tanh/i}) " "*
escapedcomputation: "{" computation "}"
NUMBER: ["-"] POSINT ["." [POSINT]] ["e"i ["-"] POSINT]

USRCONSTANT: /c[0-9]+/i
usrconstantdef: /c[0-9]+\[/i computation /]/

FILENAME: /[a-z0-9()]+/i // TODO; note that array entries can be invoked inside filenames cf. man pg. 54

%import common.WS
%ignore WS
""")


print(parser.parse(test))
