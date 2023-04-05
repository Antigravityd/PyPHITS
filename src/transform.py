from base import *

vector = Tuple(Real(), Real(), Real())
class Transform(PhitsObject): # a translation vector and a (flattened) rotation matrix representing an R^3 isometry
    syntax = {"translation": ("trans", vector, 0),
              "rotation": ("rot", Tuple(vector, vector, vector), 1),
              "rotate_first": (None, FinBij({True: 1, False: -1}), None),
              "units": (None, FinBij({"degrees": "degrees", "radians": "radians"}), None)}
    shape=((lambda: f"*TR{self.index}" if self.units == "degrees" else f"TR{self.index}",
            "trans", "rot" "rotate_first"))

#idTransform = Transform([0.0, 0.0, 0.0], [0.0 for i in range(9)])
