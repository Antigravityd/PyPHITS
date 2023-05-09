"""Regions of space that emit particles."""


import sys
from base import *
from transform import Transform
from valspec import *
from distribution import *
from cell import Cell
# TODO: global scaling factor totfact, and correlation option iscorr. Something with group_by?

# removed from projectile spec: FinBij({"all": "all"})
common = {"projectile": ("proj", List(OneOf(Particle(), Nuclide())), 0),
          "spin": (("sx", "sy", "sz"), (PosReal(), PosReal(), PosReal()), None),
          "mask": (("reg", "ntmax"), (IsA(Cell, index=True), PosInt()), None),
          "transform": ("trcl", IsA(Transform, index=True), None),
          "weight": ("wgt", PosReal(), None),
          "charge_override": ("izst", PosReal(), None),
          "counter_start": (("cnt(1)", "cnt(2)", "cnt(3)"), (PosInt(), PosInt(), PosInt()), None),
          "fissile": ("ispfs", FinBij({False: 0, "fissions": 1, "neutrons": 2}), None)
          # ibatch?
          }

semi_common = {"elevation": ("dir", OneOf(RealBetween(0.0, 1.0), FinBij({"isotropic": "all"}), IsA(AngleDistribution)), None),
               "azimuth": ("phi", PosReal(), None),
               "dispersion": ("dom", OneOf(PosReal(), FinBij({"cos^2": -1})), None),
               # "energy": ("e0", PosReal(), 1), unsupported; just use a uniform energy distribution
               "spectrum": (None, IsA(EnergyDistribution), 1)}



class Cylindrical(PhitsObject):
    name = "source"
    syntax = common | {"center": (("x0", "y0"), (Real(), Real()), None),
                        "zbounds": (("z0", "z1"), (Real(), Real()), None),
                        "radius": ("r0", PosReal(), None),
                        "cutout_radius": ("r1", PosReal(), None)} | semi_common

    shape = lambda self: ("s-type = 1", "projectile", "spin", "mask", "transform", "weight", "charge_override", "counter_start",
                          "fissile", "center", "zbounds", "radius", "cutout_radius",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))

    def restrictions(self):
        if (self.radius is None or self.radius == 0) and self.cutout_radius is not None:
            raise ValueError("Cylindrical sources that specify a cutout radius must also specify a nonzero radius;"
                             "got cutout_radius={self.cutout_radius}.")
        if self.radius is not None and self.cutout_radius is not None and self.radius < self.cutout_radius:
            raise ValueError("Cylindrical sources cannot have cutouts larger than their radius;"
                             f"got radius={self.radius} and cutout_radius={self.cutout_radius}.")

class Rectangular(PhitsObject):
    name = "source"
    syntax = common | {"xbounds": (("x0", "x1"), (Real(), Real()), None),
                       "ybounds": (("x0", "x1"), (Real(), Real()), None),
                       "zbounds": (("x0", "x1"), (Real(), Real()), None)} | semi_common

    shape = lambda self: ("s-type = 2", "projectile", "spin", "mask", "transform", "weight", "charge_override", "counter_start",
                          "fissile", "xbounds", "ybounds", "zbounds",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))



# TODO: update with correct dir definition
class Gaussian(PhitsObject):
    name = "source"
    syntax = common | {"center": (("x0", "y0", "z0"), (Real(), Real(), Real()), None),
                       "fwhms": (("x1", "y1", "z1"), (PosReal(), PosReal(), PosReal()), None)} | semi_common

    shape = lambda self: ("s-type = 3", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "center", "fwhms",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                          else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))

class GaussianPrism(PhitsObject):
    name = "source"
    syntax = common | {"center": (("x0", "y0"), (Real(), Real()), None),
                       "fwhm": ("r1", PosReal(), None),
                       "zbounds": (("z0", "z1"), (Real(), Real()), None)} | semi_common

    shape = lambda self: ("s-type = 13", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "center", "fwhm", "zbounds",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                          else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))




class Parabolic(PhitsObject):
    name = "source"

    syntax = common | {"center": (("x0", "y0"), (Real(), Real()), None),
                       "width": (("x1", "y1"), (PosReal(), PosReal()), None),
                       "zbounds": (("z0", "z1"), (Real(), Real()), None),
                       "order": ("rn", Between(2, 2147483647), None) # PHITS's default INTEGER is 32-bit; if something's bigger,
                                                                     # their ≡ 0 (mod 2) check of multiplying and dividing by 2 fails.
                       } | semi_common
    shape = lambda self: ("s-type = 7", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "center", "width", "zbounds", "order",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))

    def restrictions(self):
        if self.order is not None and self.order % 2 != 0:
            raise ValueError(f"The order of a Parabolic source must be even; got order={self.order}.") # TODO: needed?
        if self.zbounds is not None and self.zbounds[0] > self.zbounds[1]:
            raise ValueError(f"The the zbounds of a Parabolic source must be a well-formed interval; got zbounds={self.zbounds}.")

# The difference between these two in the manual is...sus
class ParabolicPrism(PhitsObject):
    name = "source"
    syntax = common | {"center": (("x0", "y0"), (Real(), Real()), None),
                       "width": ("r1", Real(), None),
                       "zbounds": (("z0", "z1"), (Real(), Real()), None),
                       "order": ("rn", Between(2, 2147483647), None)
                       } | semi_common
    shape = lambda self: ("s-type = 15", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "center", "width", "zbounds", "order",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))

    def restrictions(self):
        if self.order is not None and self.order % 2 != 0:
            raise ValueError(f"The order of a ParabolicPrism source must be even; got order={self.order}.") # TODO: needed?
        if self.zbounds is not None and self.zbounds[0] > self.zbounds[1]:
            raise ValueError(f"The the zbounds of a ParabolicPrism source must be a well-formed interval; got zbounds={self.zbounds}.")


# dir = iso not supported
class Spherical(PhitsObject):
    name = "source"
    syntax = common | {"center": (("x0", "y0", "z0"), (Real(), Real(), Real()), None),
                       "r_in": ("r1", PosReal(), None),
                       "r_out": ("r2", PosReal(), None),
                       # "elevation_bounds": (("ag1", "ag2"), (Real(), Real()), None),
                       # "azimuth_bounds": (("pg1", "pg2"), (Real(), Real()), None),
                       "elevation": ("dir", OneOf(RealBetween(0.0, 1.0), FinBij({"all": "all"}), IsA(AngleDistribution)), None),
                       # TODO: this elevation and this elevation only doesn't work if I set FinBij({"isotropic": "all"}).
                       "resample_cutoff": ("isbias", Choice10(), None),
                       "spectrum": (None, IsA(EnergyDistribution), 1)}
    shape = lambda self: ("s-type = 9", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "center", "r_in", "r_out",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "resample_cutoff", ("spectrum",))

    def restrictions(self):
        if (self.elevation == "isotropic" and self.r_in is not None and (self.r_out is None or self.r_out == 0)) \
           or (self.elevation == "isotropic" and self.r_in is not None and self.r_out is not None and self.r_in <= self.r_out):
            raise ValueError("Spherical sources with isotropic elevation must have greater inner radius than outer radius;"
                             f"got r_in={self.r_in} and r_out={self.r_out}.")


        if (self.elevation != "isotropic" and self.r_in is not None and (self.r_out is None or self.r_out == 0)) \
           or (self.elevation != "isotropic" and self.r_in is not None and self.r_out is not None and self.r_in > self.r_out):
            raise ValueError("Spherical sources that specify an inner radius must also specify a greater outer radius;"
                             f"got r_in={self.r_in} and r_out={self.r_out}.")


class Beam(PhitsObject): # I don't understand what this is trying to do
    name = "source"
    syntax = common | {"center": (("x0", "y0"), (Real(), Real()), None),
                       "eccentricity": (("x1", "y1"), (Real(), Real()), None),
                       "zbounds": (("z0", "z1"), (Real(), Real()), None),
                       "phase_gradients": (("rx", "ry"), (Real(), Real()), None),
                       "sampling": ("wem", OneOf(FinBij({"gaussian": 0}), PosReal()), None),
                       "dispersion": (("x1", "y1"), (Real(), Real()), None),
                       "angle_dispersion": (("xmrad1", "ymrad1"), (PosReal(), PosReal()), None),
                       "phase_center": (("x2", "y2"), (Real(), Real()), None),
                       "phase_angle_center": (("xmrad2", "ymrad2"), (Real(), Real()), None),
                       "positive": ("dir", FinBij({True: 1, False: -1}), None),
                       "spectrum": (None, IsA(EnergyDistribution), 1)}

    shape = ("s-type = 11", "projectile", "spin", "mask", "transform", "weight", "counter_start",
             "charge_override", "fissile", "center", "eccentricity", "zbounds", "phase_gradients", "sampling", "dispersion",
             "angle_dispersion", "phase_center", "phase_angle_center", "positive", ("spectrum",))


# decay-turtle??????


class Conical(PhitsObject):
    name = "source"
    syntax = common | {"top": (("x0", "y0", "z0"), (Real(), Real(), Real()), None),
                       "altitude": (("x1", "y1", "z1"), (Real(), Real(), Real()), None),
                       "trim": (("r0", "r1"), (Real(), Real()), None),
                       "angle": ("r2", PosReal(), None)} | semi_common
    shape = lambda self: ("s-type = 18", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "top", "altitude", "trim", "angle",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth",
                          "dispersion", ("spectrum",))



class TrianglePrism(PhitsObject):
    name = "source"
    syntax = common | {"origin": (("x0", "y0", "z0"), (Real(), Real(), Real()), None),
                       "side1": (("x1", "y1", "z1"), (Real(), Real(), Real()), None),
                       "side2": (("x2", "y2", "z2"), (Real(), Real(), Real()), None),
                       "extrusion": (("x3", "y3", "z3"), (Real(), Real(), Real()), None),
                       "attenuation": ("exa", PosReal(), None)} | semi_common
    shape = lambda self: ("s-type = 20", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "origin", "side1", "side2", "extrusion", "attenuation",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth",
                          "dispersion", ("spectrum",))



# class Grid(PhitsObject):
#     name = "source"
#     syntax = common | {"meshes": (("x-type", "y-type", "z-type"), ())}
#     required = ["projectile", "energy", "mesh"]
#     positional = ["projectile", "energy", "mesh"]
#     optional = ["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
#                 , "azimuth", "dispersion", "e0", "cutoff_behavior"]
#     ident_map = {"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
#                  "weight": "wgt", "charge_override": "izst", "fissile": "ispfs",
#                  "elevation": "dir", "azimuth": "phi", "dispersion": "dom", "energy": "e0"}
#     value_map = {"neutrons": 2, True: 1}
#     shape = ("s-type = 22", "projectile", "spin", "mask", "transform", "weight", "factor",
#              "charge_override", "fissile", "mesh", , "azimuth", "dispersion", "energy")



# class TetrahedralSource(PhitsObject): # TODO: subobjects
#     name = "source"
#     syntax = common | {"cell": ("tetreg", IsA(Cell), 2)} | semi_common
#     shape = ("s-type = 24", "projectile", "spin", "mask", "transform", "weight", "counter_start",
#              "charge_override", "fissile", "cell", "elevation", "azimuth", "dispersion", "spectrum")


# class SurfaceSource(PhitsObject):
#     name = "source"
#     syntax = common | {"surface": ("suf", IsA(Surface), 2),
#                        "cut": ("cut", IsA(Cell), 3)} | semi_common # TODO: cut sus

#     shape = ("s-type = 26", "projectile", "spin", "mask", "transform", "weight", "counter_start",
#              "charge_override", "fissile", "surface", "cut", "elevation", "azimuth", "dispersion", "spectrum")

# dump file

# user source

# class Duct(PhitsObject):
#     name = "source"
#     required = ["wall", "dl0", "dl1", "dl2", "dpf", "drd"]
#     positional = ["wall", "dl0", "dl1", "dl2", "dpf", "drd"]
#     optional = ["dxw", "dyw"]
#     def __init__(self, *args, **kwargs):

__pdoc__ = dict()
__pdoc__["builds"] = False
__pdoc__["slices"] = False
for name, cl in list(sys.modules[__name__].__dict__.items()):
    if type(cl) == type and issubclass(cl, PhitsObject) and cl != PhitsObject:
        __pdoc__[cl.__name__] = cl.__doc__ + cl.syntax_desc() if cl.__doc__ else cl.syntax_desc()
