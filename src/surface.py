import numpy as np
import sys
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
                           "transform", "P", "A", "B", "C", "D"),)

    def restrictions(self):
        if self.A == 0 and self.B == 0 and self.C == 0:
            raise ValueError("For Plane: at least one of A, B, or C must be nonzero.")


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
                           f"{self.p3[0]}", f"{self.p3[1]}", f"{self.p3[2]}"),)

    def restrictions(self):
        if self.p1[0] * (self.p2[1] - self.p3[1]) + self.p2[0] * (self.p3[1] - self.p1[1]) \
           + self.p3[0] * (self.p1[1] - self.p2[1]) == 0: # i.e. points are colinear
            raise ValueError("For PointPlane: p1, p2, and p3 must not line on a line;"
                             f" got p1={self.p1}, p2={self.p2}, and p3={self.p3}.")


class ParallelPlane(PhitsObject):
    """A plane of the form x_i = D."""
    name = "surface"
    syntax = common | {"parallel": (None, FinBij({"x": "X", "y": "Y", "z":"Z"}), 0),
                       "D": (None, Real(), 1)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", f"P{self.parallel}", "D"),)


class Sphere(PhitsObject):
    "A sphere of radius R centered on (x0, y0, z0)."
    name = "surface"
    syntax = common | {"radius": (None, PosReal(), 0),
                       "center": (None, Tuple(Real(), Real(), Real()), 1)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform",
                           "SPH", f"{self.center[0]}", f"{self.center[1]}", f"{self.center[2]}", "radius"),)


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
                           "RCC", " ".join(str(i) for i in self.center), " ".join(str(i) for i in self.height), "radius"),)

    def restrictions(self):
        if self.height == (0, 0, 0):
            raise ValueError("Cylinder must have a nonzero height vector.")

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
                           "TRC", " ".join(str(i) for i in self.center), " ".join(str(i) for i in self.height), "bottom_r", "top_r"),)

    def restrictions(self):
        if self.height == (0, 0, 0):
            raise ValueError("Cone must have a nonzero height vector.")

        if self.bottom_r <= self.top_r:
            raise ValueError("Cone must have a top radius smaller than its bottom radius;"
                             f" got bottom_r={self.bottom_r} and top_r={self.top_r}.")




class SimpleConic(PhitsObject): # ellipsoid, hyperboloid, or paraboloid parallel to an axis of the form
                   # A(x-x0)^2+B(y-y0)^2+C(z-z0)^2+2D(x-x0)+2E(y-y0)+2F(z-z0)+G = 0
    name = "surface"
    syntax = common | {"quadratic": ((None, None, None), (Real(), Real(), Real()), 0),
                       "linear": ((None, None, None), (Real(), Real(), Real()), 1),
                       "constant": (None, Real(), 2),
                       "center": ((None, None, None), (Real(), Real(), Real()), 3)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "SQ", "quadratic", "linear", "constant", "center"),)



class GeneralConic(PhitsObject): # ellipsoid, hyperboloid, or paraboloid of the form
                    # A(x-x0)^2+B(y-y0)^2+C(z-z0)^2+Dxy+Eyz+Fzx+Gx+Hy+Jz+K = 0
    name = "surface"
    syntax = common | {"quadratic": ((None, None, None), (Real(), Real(), Real()), 0),
                       "mixed": ((None, None, None), (Real(), Real(), Real()), 1),
                       "linear": ((None, None, None), (Real(), Real(), Real()), 2),
                       "constant": (None, Real(), 2)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "GQ", "quadratic", "mixed", "linear", "constant"),)

# TODO: I don't know what "skewed" means for transfomations on tori, so disabling for now.
# class Torus(PhitsObject): # torus parallel to an axis of the form
#              # (axisvar - axis0)^2/B^2 + (quadrature(<non-axis displacements>) - A)^2 - 1 = 0
#     name = "surface"
#     syntax = common | {"axis": (None, FinBij({"x": "X", "y": "Y", "z":"Z"}), 0),
#                        "center": ((None, None, None), (Real(), Real(), Real()), 1),
#                        "scales": ((None, None, None), (Real(), PosReal(), PosReal()), 2)}
#     shape = lambda self: ((f"*{self.index}" if self.reflective else
#                            (f"+{self.index}" if self.white else f"{self.index}"),
#                            f"T{self.axis}", "center", "scales"),)

#     def restrictions(self):
#         if self.scales[0] == 0 or self.scales[1] == 0 or self.scales[2] == 0:
#             raise ValueError(f"Torus's scales must be nonzero; got {self.scales}")
#         # if self.transform is not None and not self.transform.rotate_first: # skew
#         #     raise ValueError()


class Box(PhitsObject): # box formed by three vectors with tails at a given base point, or cross product of 3 intervals,
           # stored in the form x0 y0 z0 Ax Ay Az Bx By Bz Cx Cy Cz
    name = "surface"
    syntax = common | {"base": ((None, None, None), (Real(), Real(), Real()), 0),
                       "walls": (None, OrthogonalMatrix(), 1)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "BOX", "base", " ".join(" ".join(str(i) for i in j) for j in self.walls)),)


# Too complicated to generate
# class HexagonalPrism(PhitsObject):
#     name = "surface"
#     syntax = common | {"base": ((None, None, None), (Real(), Real(), Real()), 0),
#                        "height": ((None, None, None), (Real(), Real(), Real()), 1),
#                        "s1": ((None, None, None), (Real(), Real(), Real()), 2),
#                        "s2": ((None, None, None), (Real(), Real(), Real()), 3),
#                        "s3": ((None, None, None), (Real(), Real(), Real()), 4)}

#     shape = lambda self: ((f"*{self.index}" if self.reflective else
#                            (f"+{self.index}" if self.white else f"{self.index}"),
#                            "transform", "HEX", "base", "height", "s1", "s2", "s3"),)

#     def restrictions(self):
#         if self.height == (0, 0, 0) or self.s1 == (0, 0, 0) or self.s2 == (0, 0, 0) or self.s3 == (0, 0, 0):
#             raise ValueError("HexagonalPrism must have a nonzero height vector.")


class EllipticalCylinder(PhitsObject):
    name = "surface"
    syntax = common | {"center": ((None, None, None), (Real(), Real(), Real()), 0),
                       "axes": (None, OrthogonalMatrix(), 1)}
    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "REC", "center", " ".join(" ".join(str(i) for i in j) for j in self.axes)),)



class Spheroid(PhitsObject):
    name = "surface"
    syntax = common | {"focus1": ((None, None, None), (Real(), Real(), Real()), 0),
                       "focus2": ((None, None, None), (Real(), Real(), Real()), 1),
                       "major_axis": (None, Real(), 2)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "ELL", "focus1", "focus2", "major_axis"),)
    def restrictions(self):
        if self.focus1 == self.focus2:
            raise ValueError(f"Spheroid must have distinct foci; got focus1={self.focus1} and focus2={self.focus2}.")

        if self.major_axis == 0:
            raise ValueError("Spheroid must have a nonzero major axis length.")

        if self.major_axis - np.linalg.norm(np.array(self.focus1) - np.array(self.focus2)) <= 0:
            raise ValueError("Spheroid must have nonzero major axis length larger than the distance betwen its foci;"
                             f" got major_axis={self.major_axis}, focus1={self.focus1}, and focus2={self.focus2}.")


class Wedge(PhitsObject):
    name = "surface"
    syntax = common | {"tip": ((None, None, None), (Real(), Real(), Real()), 0),
                       "sides": (None, OrthogonalMatrix(), 1)}


    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "WED", "tip", " ".join(" ".join(str(i) for i in j) for j in self.sides)),)



class TetrahedronBox(PhitsObject):
    name = "surface"
    syntax = common | {"xrange": ((None, None), (Real(), Real()), 0),
                       "yrange": ((None, None), (Real(), Real()), 1),
                       "zrange": ((None, None), (Real(), Real()), 2)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                            (f"+{self.index}" if self.white else f"{self.index}"),
                            "transform", "RPP", "xrange", "yrange", "zrange"),)

    def restrictions(self):
        if self.xrange[0] >= self.xrange[1] or self.yrange[0] >= self.yrange[1] or self.zrange[0] >= self.zrange[1]:
            raise ValueError("EllipticalCylinder must have well-formed range intevals;"
                             f" got xrange={self.xrange}, yrange={self.yrange}, zrange={self.zrange}.")


surface_spec = OneOf(IsA(Plane, index=True), IsA(PointPlane, index=True), IsA(ParallelPlane, index=True),
                     IsA(Sphere, index=True), IsA(Cylinder, index=True), IsA(Cone, index=True), IsA(SimpleConic, index=True),
                     IsA(GeneralConic, index=True), IsA(Box, index=True), # IsA(Torus, index=True), IsA(HexagonalPrism, index=True),
                     IsA(EllipticalCylinder, index=True), IsA(Spheroid, index=True), IsA(Wedge, index=True), IsA(TetrahedronBox, index=True))

__pdoc__ = dict()
__pdoc__["builds"] = False
__pdoc__["slices"] = False
for name, cl in list(sys.modules[__name__].__dict__.items()):
    if type(cl) == type and issubclass(cl, PhitsObject) and cl != PhitsObject:
        __pdoc__[cl.__name__] = cl.__doc__ + cl.syntax_desc() if cl.__doc__ else cl.syntax_desc()
