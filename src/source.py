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
    mapping = common | {"center": (("x0", "y0"), (Real(), Real()), None),
                        "zbounds": (("z0", "z1"), (Real(), Real()), None),
                        "radius": ("r0", PosReal(), None),
                        "cutout_radius": ("r1", PosReal(), None)} | semi_common

    shape=("s-type = 1", "projectile", "spin", "mask", "transform", "weight", "charge_override", "counter_start",
           "fissile", "center", "zbounds", "radius", "cutout_radius", "elevation", "azimuth", "dispersion", ("spectrum",))


class Rectangular(PhitsObject):
    name = "source"
    mapping = common | {"xbounds": (("x0", "x1"), (Real(), Real()), None),
                        "ybounds": (("x0", "x1"), (Real(), Real()), None),
                        "zbounds": (("x0", "x1"), (Real(), Real()))} | semi_common

    shape=("s-type = 2", "projectile", "spin", "mask", "transform", "weight", "charge_override", "counter_start"
           "fissile", "xbounds", "ybounds", "zbounds", "elevation", "azimuth", "dispersion", ("spectrum",))



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
                          "charge_override", "fissile", "center", "fwhms", "zbounds",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                          else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))




class Parabolic(PhitsObject):
    name = "source"

    syntax = common | {"center": (("x0", "y0"), (Real(), Real()), None),
                       "width": (("x1", "y1"), (Real(), Real()), None),
                       "zbounds": (("z0", "z1"), (Real(), Real()), None),
                       "order": ("rn", Integer(), None)} | semi_common
    shape = ("s-type = 7", "projectile", "spin", "mask", "transform", "weight", "counter_start",
             "charge_override", "fissile", "center", "width", "zbounds", "order", "elevation", "azimuth",
             "dispersion", ("spectrum",))

# The difference between these two in the manual is...sus
class ParabolicPrism(PhitsObject):
    name = "source"
    syntax = common | {"center": (("x0", "y0"), (Real(), Real()), None),
                       "width": ("r1", Real(), None),
                       "zbounds": (("z0", "z1"), (Real(), Real()), None),
                       "order": ("rn", Integer(), None)} | semi_common
    shape = ("s-type = 15", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                  "charge_override", "fissile", "center", "width", "zbounds", "order", "elevation", "azimuth",
                  "dispersion", ("spectrum",))


# dir = iso not supported
class Spherical(PhitsObject):
    name = "source"
    syntax = common | {"center": (("x0", "y0", "z0"), (Real(), Real(), Real()), None),
                       "r_in": ("r1", Real(), None),
                       "r_out": ("r2", Real(), None),
                       # "elevation_bounds": (("ag1", "ag2"), (Real(), Real()), None),
                       # "azimuth_bounds": (("pg1", "pg2"), (Real(), Real()), None),
                       "elevation": ("dir", OneOf(PosReal(), FinBij({"isotropic": "all"}), IsA(AngleDistribution)), None),
                       "resample_cutoff": ("isbias", Choice10(), None),
                       "spectrum": (None, IsA(EnergyDistribution), 1)}
    shape = ("s-type = 9", "projectile", "spin", "mask", "transform", "weight", "counter_start",
             "charge_override", "fissile", "center", "r_in", "r_out", "elevation", "resample_cutoff", ("spectrum",))



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
                       "positive": ("dir", Choice10(true=1, false=-1), None)}

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
    shape = ("s-type = 18", "projectile", "spin", "mask", "transform", "weight", "counter_start",
             "charge_override", "fissile", "top", "altitude", "trim", "angle", "elevation", "azimuth",
             "dispersion", ("spectrum",))



class TrianglePrism(PhitsObject):
    name = "source"
    syntax = common | {"origin": (("x0", "y0", "z0"), (Real(), Real(), Real()), None),
                       "side1": (("x1", "y1", "z1"), (Real(), Real(), Real()), None),
                       "side2": (("x2", "y2", "z2"), (Real(), Real(), Real()), None),
                       "extrusion": (("x3", "y3", "z3"), (Real(), Real(), Real()), None),
                       "attenuation": ("exa", PosReal(), None)} | semi_common
    shape = ("s-type = 20", "projectile", "spin", "mask", "transform", "weight", "counter_start",
             "charge_override", "fissile", "origin", "side1", "side2", "extrusion", "attenuation", "elevation", "azimuth",
             "dispersion", ("spectrum",))



# class Grid(PhitsObject):
#     name = "source"
#     syntax = common | {"meshes": (("x-type", "y-type", "z-type"), ())}
#     required = ["projectile", "energy", "mesh"]
#     positional = ["projectile", "energy", "mesh"]
#     optional = ["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
#                 "elevation", "azimuth", "dispersion", "e0", "cutoff_behavior"]
#     ident_map = {"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
#                  "weight": "wgt", "charge_override": "izst", "fissile": "ispfs",
#                  "elevation": "dir", "azimuth": "phi", "dispersion": "dom", "energy": "e0"}
#     value_map = {"neutrons": 2, True: 1}
#     shape = ("s-type = 22", "projectile", "spin", "mask", "transform", "weight", "factor",
#              "charge_override", "fissile", "mesh", "elevation", "azimuth", "dispersion", "energy")



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
