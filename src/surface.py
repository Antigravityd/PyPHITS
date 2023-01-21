import numpy as np
from base import *
from transform import *



class Plane(PhitsObject): # Planar surface stored as Ax+By+Cz-D = 0.
    def __init__(self, *args, **kwargs):
        super().no_hash.add("inside")
        if len(args) == 4:
            self.name = "surface"
            self.required = ["A", "B", "C", "D"]
            self.positional = ["A", "B", "C", "D"]
            self.optional = ["reflective", "white", "transform", "inside"]
            self.shape = ((lambda self: f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "P", "A", "B", "C", "D"))
            super().__init__(*args, **kwargs)
        elif len(args) == 3:
            self.name = "surface"
            self.required = ["p1", "p2", "p3"]
            self.positional = ["p1", "p2", "p3"]
            self.optional = ["reflective", "white", "transform", "inside"]
            self.shape = ((lambda self: f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"), "transform", "P", "p1", "p2", "p3"),)
            super().__init__(*args, **kwargs)
        else:
            self.name = "surface"
            self.required = ["parallel", "D"]
            self.positional = ["parallel", "D"]
            self.optional = ["reflective", "white", "transform", "inside"]
            self.shape = ((lambda self: f"*{self.index}" if self.reflective else
                      (f"+{self.index}" if self.white else f"{self.index}"),
                      "transform",
                      lambda self: f"P{self.parallel}", "D"))
            super().__init__(*args, **kwargs)


class Sphere(PhitsObject): # Spherical surface stored as (x-x0)^2+(y-y0)^2+(z-z0)^2-R^2=0
    def __init__(self, *args, **kwargs):
        super().no_hash.add("inside")
        if len(args) == 1:
            self.name = "surface"
            self.required = ["R"]
            self.positional = ["R"]
            self.optional = ["reflective", "white", "transform", "inside"]
            self.shape = ((lambda self: f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"), "transform", "SO", "R"),)
            super().__init__(*args, **kwargs)
        elif len(args) == 2:
            self.name = "surface"
            self.required = ["center", "R"]
            self.positional = ["center", "R"]
            self.optional = ["reflective", "white", "transform", "inside"]
            self.shape = ((lambda self: f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"), "transform", "S", "center", "R"),)
            super().__init__(*args, **kwargs)
        elif len(args) == 3:
            self.name = "surface"
            self.required = ["coordinate", "on", "R"]
            self.positional = ["coordinate", "on" "R"]
            self.optional = ["reflective", "white", "transform", "inside"]
            self.shape = ((lambda self: f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", lambda self: f"S{self.on}", "coordinate" "R"),)
            super().__init__(*args, **kwargs)


class Cylinder(PhitsObject):
    def __init__(self, *args, **kwargs):
        super().no_hash.add("inside")
        symb = args[0]
        if len(args) == 2:
            self.name = "surface"
            self.required = ["on", "R"]
            self.positional = ["on", "R"]
            self.optional = ["reflective", "white", "transform", "inside"]
            self.shape = ((lambda self: f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", lambda self: f"C{self.on}", "R"),)
            super().__init__(*args, **kwargs)
        elif len(args) == 3:
            self.name = "surface"
            self.required = ["parallel", "center", "R"]
            self.positional = ["parallel", "center", "R"]
            self.optional = ["reflective", "white", "transform", "inside"]
            self.shape = ((lambda self: f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", lambda self: f"C/{self.parallel}", "center", "R"),)
            super().__init__(*args, **kwargs)


class Cone(PhitsObject):
    def __init__(self, *args, **kwargs):
        super().no_hash.add("inside")
        if isinstance(args[1], list):
            self.name = "surface"
            self.required = ["on", "coordinate", "t_squared", "sheet"]
            self.positional = ["on","coordinate", "t_squared", "sheet"]
            self.optional = ["reflective", "white", "transform", "inside"]
            self.shape = ((lambda self: f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", lambda self: f"K{self.on}", "coordinate", "t_squared", "sheet"),)
            self.value_map = {"upper": 1.0, "lower": -1.0}
            super().__init__(*args, **kwargs)

        else:
            self.name = "surface"
            self.required = ["parallel", "coordinate", "t_squared", "sheet"]
            self.positional = ["parallel","coordinate", "t_squared", "sheet"]
            self.optional = ["reflective", "white", "transform", "inside"]
            self.shape = ((lambda self: f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", lambda self: f"K/{self.parallel}", "coordinate", "t_squared", "sheet"),),
            self.value_map = {"upper": 1.0, "lower": -1.0}
            super().__init__(*args, **kwargs)


class SimpleConic(PhitsObject): # ellipsoid, hyperboloid, or paraboloid parallel to an axis of the form
                   # A(x-x0)^2+B(y-y0)^2+C(z-z0)^2+2D(x-x0)+2E(y-y0)+2F(z-z0)+G = 0
    name = "surface"
    required = ["A", "B", "C", "D", "E", "F", "G", "center"]
    positional = ["A", "B", "C", "D", "E", "F", "G", "center"]
    optional = ["reflective", "white", "transform", "inside"]
    shape = ((lambda self: f"*{self.index}" if self.reflective else
              (f"+{self.index}" if self.white else f"{self.index}"),
              "transform", "SQ", "A", "B", "C", "D", "E", "F", "G", "center"),)
    def __init__(self, *args, **kwargs):
        super().no_hash.add("inside")



class GeneralConic(PhitsObject): # ellipsoid, hyperboloid, or paraboloid of the form
                    # A(x-x0)^2+B(y-y0)^2+C(z-z0)^2+Dxy+Eyz+Fzx+Gx+Hy+Jz+K = 0
    name = "surface"
    required = ["A", "B", "C", "D", "E", "F", "G", "H", "J", "K"]
    positional = ["A", "B", "C", "D", "E", "F", "G", "H", "J", "K"]
    optional = ["reflective", "white", "transform", "inside"]
    shape = ((lambda self: f"*{self.index}" if self.reflective else
              (f"+{self.index}" if self.white else f"{self.index}"),
              "transform", "GQ", "A", "B", "C", "D", "E", "F", "G", "H", "J", "K"),)
    def __init__(self, *args, **kwargs):
        super().no_hash.add("inside")



class Torus(PhitsObject): # torus parallel to an axis of the form
             # (axisvar - axis0)^2/B^2 + (quadrature(<non-axis displacements>) - A)^2 - 1 = 0
    name = "surface"
    required = ["parallel", "center", "A", "B", "C"]
    positional = ["parallel", "center", "A", "B", "C"]
    optional = ["reflective", "white", "transform", "inside"]
    shape = ((lambda self: f"*{self.index}" if self.reflective else
              (f"+{self.index}" if self.white else f"{self.index}"),
              "transform", lambda self: f"T{self.parallel}", "center", "A", "B", "C"),)
    def __init__(self, *args, **kwargs):
        super().no_hash.add("inside")



class Box(PhitsObject): # box formed by three vectors with tails at a given base point, or cross product of 3 intervals,
           # stored in the form x0 y0 z0 Ax Ay Az Bx By Bz Cx Cy Cz
    name = "surface"
    required = ["base", "s1", "s2", "s3"]
    positional = ["base", "s1", "s2", "s3"]
    optional = ["reflective", "white", "transform", "inside"]
    shape = ((lambda self: f"*{self.index}" if self.reflective else
              (f"+{self.index}" if self.white else f"{self.index}"),
              "transform", "BOX", "base", "s1", "s2", "s3"),)
    def __init__(self, *args, **kwargs):
        super().no_hash.add("inside")


class HexagonalPrism(PhitsObject):
    name = "surface"
    required = ["base", "height", "s1", "s2", "s3"]
    positional = ["base", "height", "s1", "s2", "s3"]
    optional = ["reflective", "white", "transform", "inside"]
    shape = ((lambda self: f"*{self.index}" if self.reflective else
              (f"+{self.index}" if self.white else f"{self.index}"),
              "transform", "HEX", "base", "height", "s1", "s2", "s3"))
    def __init__(self, *args, **kwargs):
        super().no_hash.add("inside")


class ElipticalCylinder(PhitsObject):
    name = "surface"
    required = ["center", "height", "major_axis", "minor_axis"]
    positional = ["center", "height", "major_axis", "minor_axis"]
    optional = ["reflective", "white", "transform", "inside"]
    shape = ((lambda self: f"*{self.index}" if self.reflective else
              (f"+{self.index}" if self.white else f"{self.index}"),
              "transform", "REC", "center", "height", "major_axis", "minor_axis"),)
    def __init__(self, *args, **kwargs):
        super().no_hash.add("inside")



class TruncatedCone(PhitsObject):
    name = "surface"
    required = ["center", "height", "r_top", "r_bottom"]
    positional = ["center", "height", "r_top", "r_bottom"]
    optional = ["reflective", "white", "transform", "inside"]
    shape = ((lambda self: f"*{self.index}" if self.reflective else
              (f"+{self.index}" if self.white else f"{self.index}"),
              "transform", "TRC", "center", "height", "r_top", "r_bottom"))
    def __init__(self, *args, **kwargs):
        super().no_hash.add("inside")



class Spheroid(PhitsObject):
    name = "surface"
    required = ["focus_1", "focus_2" "major_axis"]
    positional = ["focus_1", "focus_2" "major_axis"]
    optional = ["reflective", "white", "transform", "inside"]
    shape = ((lambda self: f"*{self.index}" if self.reflective else
              (f"+{self.index}" if self.white else f"{self.index}"),
              "transform", "ELL", "focus_1", "focus_2", "major_axis"),)
    def __init__(self, *args, **kwargs):
        super().no_hash.add("inside")


class Wedge(PhitsObject):
    name = "surface"
    required = ["tip", "s1", "s2", "height"]
    positional = ["tip", "s1", "s2", "height"]
    optional = ["reflective", "white", "transform", "inside"]
    shape = ((lambda self: f"*{self.index}" if self.reflective else
              (f"+{self.index}" if self.white else f"{self.index}"),
              "transform", "WED", "tip", "s1", "s2", "height"),)
    def __init__(self, *args, **kwargs):
        super().no_hash.add("inside")
