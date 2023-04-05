import numpy as np
from base import *
from transform import *


common = {"reflective": (None, Choice10(), None),
          "white": (None, Choice10(), None),
          "transform": (None, IsA(Transform, index=True), None),
          "inside": (None, Choice10(), None)}


class Plane(PhitsObject):
    """A plane of the form Ax + By + Cz - D = 0."""
    name = "surface"
    syntax = common | {"A": (None, Real(), 0),
                       "B": (None, Real(), 1),
                       "C": (None, Real(), 2),
                       "D": (None, Real(), 3)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "P", "A", "B", "C", "D"))


# TODO: consider obliterating the next 2
class PointPlane(PhitsObject):
    """A plane specified by three points."""
    name = "surface"
    syntax = common | {"p1": (None, Tuple(Real(), Real(), Real()), 0),
                       "p2": (None, Tuple(Real(), Real(), Real()), 1),
                       "p3": (None, Tuple(Real(), Real(), Real()), 2)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "P",
                           f"{self.p1[0]}", f"{self.p1[1]}", f"{self.p1[2]}",
                           f"{self.p2[0]}", f"{self.p2[1]}", f"{self.p2[2]}",
                           f"{self.p3[0]}", f"{self.p3[1]}", f"{self.p3[2]}"))

class ParallelPlane(PhitsObject):
    """A plane of the form x_i = D."""
    name = "surface"
    syntax = common | {"parallel": (None, FinBij({"x": "X", "y": "Y", "z":"Z"}), 0),
                       "D": (None, Real(), 1)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", f"P{self.parallel}", "D"))


class Sphere(PhitsObject):
    "A sphere of radius R centered on (x0, y0, z0)."
    name = "surface"
    syntax = common | {"radius": (None, PosReal(), 0),
                       "center": (None, Tuple(Real(), Real(), Real()), 1)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform",
                           "SPH", f"{self.center[0]}", f"{self.center[1]}", f"{self.center[2]}" "R"))


class Cylinder(PhitsObject):
    """A right-circular cylinder with center of the bottom face (x_0, y_0, z_0), height vector from the bottom to top face (H_x, H_y, H_z),
    and radius R."""
    name = "surface"
    syntax = common | {"center": (None, Tuple(Real(), Real(), Real()), 0),
                       "height": (None, Tuple(Real(), Real(), Real()), 1),
                       "radius": (None, PosReal(), 2)}
    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform",
                           "RCC", " ".join(self.center), " ".join(self.height), "radius"))

class Cone(PhitsObject):
    """A truncated right-angle cone with bottom-face center (x_0, y_0, z_0), height vector (H_x, H_y, H_z), and bottom and top radii
    R_1 and R_2 respectively."""
    name = "surface"
    syntax = common | {"center": (None, Tuple(Real(), Real(), Real()), 0),
                       "height": (None, Tuple(Real(), Real(), Real()), 1),
                       "bottom_r": (None, PosReal(), 2),
                       "top_r": (None, PosReal(), 3)}
    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform",
                           "TRC", " ".join(self.center), " ".join(self.height), "bottom_r", "top_r"))


class SimpleConic(PhitsObject): # ellipsoid, hyperboloid, or paraboloid parallel to an axis of the form
                   # A(x-x0)^2+B(y-y0)^2+C(z-z0)^2+2D(x-x0)+2E(y-y0)+2F(z-z0)+G = 0
    name = "surface"
    syntax = common | {"quadratic": ((None, None, None), (Real(), Real(), Real()), 0),
                       "linear": ((None, None, None), (Real(), Real(), Real()), 1),
                       "constant": (None, Real(), 2),
                       "center": ((None, None, None), (Real(), Real(), Real()), 3)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "SQ", "quadratic", "linear", "constant", "center"))



class GeneralConic(PhitsObject): # ellipsoid, hyperboloid, or paraboloid of the form
                    # A(x-x0)^2+B(y-y0)^2+C(z-z0)^2+Dxy+Eyz+Fzx+Gx+Hy+Jz+K = 0
    name = "surface"
    syntax = common | {"quadratic": ((None, None, None), (Real(), Real(), Real()), 0),
                       "mixed": ((None, None, None), (Real(), Real(), Real()), 1)
                       "linear": ((None, None, None), (Real(), Real(), Real()), 2),
                       "constant": (None, Real(), 2)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "GQ", "quadratic", "mixed", "linear", "constant"))


class Torus(PhitsObject): # torus parallel to an axis of the form
             # (axisvar - axis0)^2/B^2 + (quadrature(<non-axis displacements>) - A)^2 - 1 = 0
    name = "surface"
    syntax = common | {"axis": (None, FinBij({"x": "X", "y": "Y", "z":"Z"}), 0),
                       "center": ((None, None, None), (Real(), Real(), Real()), 1),
                       "scales": ((None, None, None), (Real(), Real(), Real()), 1)}
    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", f"T{self.axis}", "center", "scales"))



class Box(PhitsObject): # box formed by three vectors with tails at a given base point, or cross product of 3 intervals,
           # stored in the form x0 y0 z0 Ax Ay Az Bx By Bz Cx Cy Cz
    name = "surface"
    syntax = common | {"base": ((None, None, None), (Real(), Real(), Real()), 0),
                       "s1": ((None, None, None), (Real(), Real(), Real()), 1),
                       "s2": ((None, None, None), (Real(), Real(), Real()), 2),
                       "s3": ((None, None, None), (Real(), Real(), Real()), 3)}
    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "BOX", "base", "s1", "s2", "s3"),)


class HexagonalPrism(PhitsObject):
    name = "surface"
    syntax = common | {"base": ((None, None, None), (Real(), Real(), Real()), 0),
                       "height": ((None, None, None), (Real(), Real(), Real()), 1),
                       "s1": ((None, None, None), (Real(), Real(), Real()), 2),
                       "s2": ((None, None, None), (Real(), Real(), Real()), 3),
                       "s3": ((None, None, None), (Real(), Real(), Real()), 4)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "HEX", "base", "height", "s1", "s2", "s3"))


class ElipticalCylinder(PhitsObject):
    name = "surface"
    syntax = common | {"center": ((None, None, None), (Real(), Real(), Real()), 0),
                       "height": ((None, None, None), (Real(), Real(), Real()), 1),
                       "major_axis": ((None, None, None), (Real(), Real(), Real()), 2),
                       "minor_axis": ((None, None, None), (Real(), Real(), Real()), 3)}
    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "REC", "center", "height", "major_axis", "minor_axis"),)




class Spheroid(PhitsObject):
    name = "surface"
    syntax = common | {"focus1": ((None, None, None), (Real(), Real(), Real()), 0),
                       "focus2": ((None, None, None), (Real(), Real(), Real()), 1),
                       "major_axis": (None, Real(), 2)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "ELL", "focus1", "focus2", "major_axis"),)


class Wedge(PhitsObject):
    name = "surface"
    syntax = common | {"tip": ((None, None, None), (Real(), Real(), Real()), 0),
                       "s1": ((None, None, None), (Real(), Real(), Real()), 1),
                       "s2": ((None, None, None), (Real(), Real(), Real()), 2),
                       "height": ((None, None, None), (Real(), Real(), Real()), 3)}

    shape = ((lambda self: f"*{self.index}" if self.reflective else
              (f"+{self.index}" if self.white else f"{self.index}"),
              "transform", "WED", "tip", "s1", "s2", "height"),)


class TetrahedronBox(PhitsObject):
    name = "surface"
    syntax = common | {"xrange": ((None, None), (Real(), Real()), 0),
                       "yrange": ((None, None), (Real(), Real()), 1),
                       "zrange": ((None, None), (Real(), Real()), 2)}

    shape = ((lambda self: f"*{self.index}" if self.reflective else
              (f"+{self.index}" if self.white else f"{self.index}"),
              "transform", "RPP", "xrange", "yrange", "zrange"),)
