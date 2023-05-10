from base import *
from numpy.linalg import det
import numpy as np
import sys

vector = Tuple(Real(), Real(), Real())
class Transform(PhitsObject): #
    """An \\(\\mathbb{R}^3\\) isometry represented by a translation vector and a rotation matrix.
    """
    name = "transform"
    syntax = {"translation": (None, vector, 0),
              "rotation": (None, vector, 1),
              "rotate_first": (None, FinBij({True: 2, False: -2}), None, -2),
              "units": (None, FinBij({"degrees": "degrees", "radians": "radians"}), None)}
    shape = lambda self: ((f"*TR{self.index}" if self.units == "degrees" else f"TR{self.index}",
                           " ".join(str(i) for i in self.translation),
                           " ".join(str(i) for i in self.rotation),
                           "0 0 0 0 0 0",
                           "rotate_first"),)




#idTransform = Transform([0.0, 0.0, 0.0], [0.0 for i in range(9)])

__pdoc__ = dict()
__pdoc__["builds"] = False
__pdoc__["slices"] = False
for name, cl in list(sys.modules[__name__].__dict__.items()):
    if type(cl) == type and issubclass(cl, PhitsObject) and cl != PhitsObject:
        __pdoc__[cl.__name__] = cl.__doc__ + cl.syntax_desc() if cl.__doc__ else cl.syntax_desc()
