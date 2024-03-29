from math import *
import numpy as np
from pyphits import *

# A dead whale, modelled by concentric spherical annuli of an increasingly less-concentrated carbon-water mixture,
# has swallowed an orphaned 241Am RTG source and is sinking into the ocean.
# Find how much the water will heat, to guide recovery efforts.

# Construct objects via list comprehensions
mats = [Material([("C", i/100), ("H", (2/3)*(100-i)/100), ("O", (1/3)*(100-i)/100)]) for i in range(55,100)]

# Or loops
cells = [Cell(regions=[Sphere(1, center=(0, 0, 0))], material=mats[0], -1)]
for i in range(1,11):
    new = Cell([Sphere(i+1), Sphere(i, inside=False)], mats[i], -1)
    cells.append(new)

ocean = Void([Sphere(11, inside=False)])

escape_e = DumpFluence(cells[-1], ocean, 4*pi/3, [1, 8], "flux")


cells.append(ocean)

source = Cylindrical("241Am", 2.2, fissile="neutrons", bounds=(-0.25,0.25), r_out=0.3)


# Capture input for further analysis in Python as your choice of many common data formats,
# or render and return the .eps
# Pass raw input data in case of an error or unimplemented feature
inp = run_phits(cells, source, escape_e, raw="$ This could be e.g. an [Elastic Option] section text.\n",
                parameters={"negs": 1, "e-mode": 2})
print(inp)
