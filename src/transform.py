from base import *
class Transform(PhitsObject): # a translation vector and a (flattened) rotation matrix representing an R^3 isometry
    def __init__(self, *args, **kwargs):
        super().__init__("transform", required=["trans", "rot"], positional=["trans", "rot"],
                         optional=["rotate_first", "units"],
                         shape=((lambda: f"*TR{self.index}" if self.units == "degrees" else f"TR{self.index}",
                                 "trans", "rot" "rotate_first")),
                         value_map={True: 1, False: -1}, *args, **kwargs)

#idTransform = Transform([0.0, 0.0, 0.0], [0.0 for i in range(9)])
