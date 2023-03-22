# Sources will be specified by an iterable whose elements are a tuple  (Source(), weight) where the last two are optional.
# totfact and iscorr may be specified in the general parameters section.

# Currently, no a-type mesh support for the elevation angles.
from base import *
from cell import Cell
from transform import Transform

# Sources will be specified by an iterable whose elements are a tuple  (Source(), weight) where the last two are optional.
# totfact and iscorr may be specified in the general parameters section.

# Currently, no a-type mesh support for the elevation angles.
from base import *
from collections import namedtuple

common = {"projectile": ("proj", FinBij({}), 0),
          "spin": (("sx", "sy", "sz"), (PosReal(), PosReal(), PosReal()), None),
          "mask": (("reg", "ntmax"), (IsA(Cell, index=True), PosInt()), None),
          "transform": ("trcl", IsA(Transform, index=True), None),
          "weight": ("wgt", PosReal(), None),
          "charge_override": ("izst", PosReal(), None),
          "counter_start": (("cnt(1)", "cnt(2)", "cnt(3)"), (PosInt(), PosInt(), PosInt()), None),
          "fissile": ("ispfs", FinBij({False: 0, "fissions": 1, "neutrons": 2}), None)
          # ibatch?
          }

semi_common = {"elevation": ("dir", OneOf(PosReal(), FinBij({"isotropic": "all"}), IsA(AngleDistribution)), None) #TODO:isa, oneof
               "azimuth": ("phi", PosReal(), None),
               "dispersion": ("dom", OneOf(PosReal(), FinBij({"cos^2": -1})), None),
               "energy": ("e0", PosReal(), 1), # TODO: make this work
               "spectrum": ("e-type", IsA(EnergyDistribution), 1)}



class Cylindrical(PhitsObject):
    mapping = common | {"center": (("x0", "y0"), (0.0, 0.0), (posreal, posreal)),
                        "zbounds": (("z0", "z1"), (0.0, 0.0), (posreal, posreal)),
                        "radius": ("r0", 0.0, posreal),
                        "cutout_radius": ("r1", 0.0, posreal)} | semi_common


class Rectangular(PhitsObject):
    mapping = common | {"xbounds": (("x0", "x1"), (0.0, 0.0), (posreal, posreal)),
                        "ybounds": (("x0", "x1"), (0.0, 0.0), (posreal, posreal)),
                        "zbounds": (("x0", "x1"), (0.0, 0.0), (posreal, posreal))} | semi_common





class Cylindrical(PhitsObject):
    name = "source"
    required=["projectile", "energy"]
    positional=["projectile", "energy"]
    optional=["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
              "center", "bounds", "r_out", "r_in", "elevation", "azimuth", "dispersion"]
    ident_map={"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
               "weight": "wgt", "charge_override": "izst", "fissile": "ispfs", "center": ("x0", "y0"),
               "bounds": ("z0", "z1"), "r_out": "r0", "r_in": "r1", "elevation": "dir", "azimuth": "phi",
               "dispersion": "dom", "energy": "e0", "projectile": "proj"}
    value_map={"neutrons": 2, True: 1}
    shape=("s-type = 1", "projectile", "spin", "mask", "transform", "weight", "factor", "charge_override",
           "fissile", "center", "bounds", "r_out", "r_in", "elevation", "azimuth", "dispersion", "energy")



class Rectangular(PhitsObject):
    name = "source"
    required=["projectile", "energy"]
    positional=["projectile", "energy"]
    optional=["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
              "xbounds", "ybounds", "zbounds", "elevation", "azimuth", "dispersion"]
    ident_map={"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
               "weight": "wgt", "charge_override": "izst", "fissile": "ispfs", "xbounds": ("x0", "x1"),
               "ybounds": ("y0", "y1"), "zbounds": ("z0", "z1"), "elevation": "dir", "azimuth": "phi",
               "dispersion": "dom", "energy": "e0"}
    value_map={"neutrons": 2, True: 1}
    shape=("s_type = 2", "projectile", "spin", "mask", "transform", "weight", "factor", "charge_override",
           "fissile", "xbounds", "ybounds", "zbounds", "elevation", "azimuth", "dispersion", "energy")



class Gaussian(PhitsObject):
    def __init__(self, *args, **kwargs):
        if "zbounds" in kwargs:
            self.name = "source"
            self.required = ["projectile", "energy"]
            self.positional = ["projectile", "energy"]
            self.optional = ["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
                      "center", "fwhms", "zbounds", "elevation", "azimuth", "dispersion"]
            self.ident_map = {"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
                       "weight": "wgt", "charge_override": "izst", "fissile": "ispfs",
                       "center": ("x0", "y0"), "fwhms": "r1", "zbounds": ("z0", "z1"),
                       "elevation": "dir", "azimuth": "phi", "dispersion": "dom", "energy": "e0"}
            self.value_map = {"neutrons": 2, True: 1}
            self.shape = ("s_type = 13", "projectile", "spin", "mask", "transform", "weight", "factor",
                   "charge_override", "fissile", "center", "fwhms", "zbounds", "elevation", "azimuth",
                   "dispersion", "energy")
            super().__init__(*args, **kwargs)
        else:
            self.name = "source"
            self.required = ["projectile", "energy"]
            self.positional = ["projectile", "energy"]
            self.optional = ["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
                      "center", "fwhms", "zbounds", "elevation", "azimuth", "dispersion"]
            self.ident_map = {"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
                              "weight": "wgt", "charge_override": "izst", "fissile": "ispfs",
                              "center": ("x0", "y0", "z0"), "fwhms": ("x1", "y1", "z1"),
                              "elevation": "dir", "azimuth": "phi", "dispersion": "dom", "energy": "e0"}
            self.value_map = {"neutrons": 2, True: 1}
            self.shape = ("s_type = 3", "projectile", "spin", "mask", "transform", "weight", "factor",
                          "charge_override", "fissile", "center", "fwhms", "elevation", "azimuth",
                          "dispersion", "energy")
            super().__init__(*args, **kwargs)




class Parabolic(PhitsObject):
    # Thanks to lazy typing, xyz or x-y is determined by dimension of center
    def __init__(self, *args, **kwargs):
        if width in kwargs and len(kwargs["width"]) == 2:
            self.name = "source"
            self.required = ["projectile", "energy"]
            self.positional = ["projectile", "energy"]
            self.optional = ["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
                             "center", "width", "zbounds", "order", "elevation", "azimuth", "dispersion"]
            self.ident_map = {"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
                              "weight": "wgt", "charge_override": "izst", "fissile": "ispfs",
                              "center": ("x0", "y0"), "width": "r1", "zbounds": ("z0", "z1"), "order": "rn",
                              "elevation": "dir", "azimuth": "phi", "dispersion": "dom", "energy": "e0"}
            self.value_map = {"neutrons": 2, True: 1}
            self.shape = ("s_type = 15", "projectile", "spin", "mask", "transform", "weight", "factor",
                          "charge_override", "fissile", "center", "width", "zbounds", "order", "elevation", "azimuth",
                          "dispersion", "energy")
            super().__init__(*args, **kwargs)
        else:
            self.name = "source"
            self.required = ["projectile", "energy"]
            self.positional = ["projectile", "energy"]
            self.optional = ["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
                      "center", "width", "zbounds", "order", "elevation", "azimuth", "dispersion"]
            self.ident_map = {"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
                       "weight": "wgt", "charge_override": "izst", "fissile": "ispfs",
                       "center": ("x0", "y0", "z0"), "width": ("x1", "y1", "z1"), "order": "rn",
                       "elevation": "dir", "azimuth": "phi", "dispersion": "dom", "energy": "e0"}
            self.value_map = {"neutrons": 2, True: 1}
            self.shape = ("s_type = 7", "projectile", "spin", "mask", "transform", "weight", "factor",
                   "charge_override", "fissile", "center", "width", "zbounds", "order", "elevation", "azimuth",
                   "dispersion", "energy")
            super().__init__(*args, **kwargs)



class Spherical(PhitsObject):
    name = "source"
    required = ["projectile", "energy"]
    positional = ["projectile", "energy"]
    optional = ["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
              "center", "r_in", "r_out", "direction", "elevation_bounds", "azimuth_bounds",
              "cutoff_behavior"]
    ident_map = {"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
               "weight": "wgt", "charge_override": "izst", "fissile": "ispfs",
               "center": ("x0", "y0", "z0"), "r_in": "r1", "r_out": "r2", "direction": "dir",
               "elevation_bounds": ("ag1", "ag2"), "azimuth_bounds": ("pg1", "pg2"),
               "cutoff_behavior": "isbias", "energy": "e0"}
    value_map = {"neutrons": 2, True: 1, "ignore": 0, "resample": 1, "outward": 1.0, "inward": -1.0,
               "isotropic": "all", "iso": "all", "cosine": "-all"}
    shape = ("s_type = 9", "projectile", "spin", "mask", "transform", "weight", "factor",
           "charge_override", "fissile", "center", "r_in", "r_out", "direction", "elevation_bounds",
           "azimuth_bound", "cutoff_behavior", "energy")


class Beam(PhitsObject): # I don't understand what this is trying to do
    name = "source"
    required = ["projectile", "energy"]
    positional = ["projectile", "energy"]
    optional = ["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
              "center", "zbounds", "gradients", "emittance", "widths", "stdevs",
              "phase_center", "phase_angle_center", "sign"]
    ident_map = {"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
               "weight": "wgt", "charge_override": "izst", "fissile": "ispfs",
               "center": ("x0", "y0"), "zbounds": ("z0", "z1"), "gradients": ("rx", "ry"),
               "emittance": "wem", "widths": ("x1", "y1"), "stdevs": ("xmrad1", "ymrad1"),
               "phase_center": ("x2", "y2"), "phase_angle_maxes": ("xmrad2", "ymrad2"), "sign": "dir", "energy": "e0"}
    value_map = {"neutrons": 2, True: 1}
    shape = ("s_type = 11", "projectile", "spin", "mask", "transform", "weight", "factor",
           "charge_override", "fissile", "center", "zbounds", "gradients", "emittance", "widths",
           "stdevs", "phase_center", "phase_angle_center", "sign")


class Conical(PhitsObject):
    name = "source"
    required = ["projectile", "energy"]
    positional = ["projectile", "energy"]
    optional = ["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
                "top", "altitude", "trim", "elevation", "azimuth", "dispersion"]
    ident_map = {"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
                 "weight": "wgt", "charge_override": "izst", "fissile": "ispfs",
                 "top": ("x0", "y0", "z0"), "altitude": ("x1", "y1", "z1"), "trim": ("r0", "r1"),
                 "elevation": "dir", "azimuth": "phi", "dispersion": "dom", "energy": "e0"}
    value_map = {"neutrons": 2, True: 1}
    shape = ("s_type = 18", "projectile", "spin", "mask", "transform", "weight", "factor",
           "charge_override", "fissile", "top", "altitude", "trim", "elevation", "azimuth",
           "dispersion", "energy")



class Prism(PhitsObject):
    name = "source"
    required = ["projectile", "energy"]
    positional = ["projectile", "energy"]
    optional = ["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
                "origin", "side1", "side2", "extrusion", "attenuation", "elevation", "azimuth", "dispersion"]
    ident_map = {"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
                 "weight": "wgt", "charge_override": "izst", "fissile": "ispfs",
                 "origin": ("x0", "y0", "z0"), "side1": ("x1", "y1", "z1"), "side2": ("x2", "y2", "z2"),
                 "extrusion": ("x3", "y3", "z3"), "attenuation": "exa",
                 "elevation": "dir", "azimuth": "phi", "dispersion": "dom", "energy": "e0"}
    value_map = {"neutrons": 2, True: 1}
    shape = ("s_type = 20", "projectile", "spin", "mask", "transform", "weight", "factor",
             "charge_override", "fissile", "origin", "side1", "side2", "extrusion", "attenuation", "elevation", "azimuth",
             "dispersion", "energy")



class Grid(PhitsObject):
    name = "source"
    required = ["projectile", "energy", "mesh"]
    positional = ["projectile", "energy", "mesh"]
    optional = ["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
                "elevation", "azimuth", "dispersion", "e0", "cutoff_behavior"]
    ident_map = {"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
                 "weight": "wgt", "charge_override": "izst", "fissile": "ispfs",
                 "elevation": "dir", "azimuth": "phi", "dispersion": "dom", "energy": "e0"}
    value_map = {"neutrons": 2, True: 1}
    shape = ("s_type = 22", "projectile", "spin", "mask", "transform", "weight", "factor",
             "charge_override", "fissile", "mesh", "elevation", "azimuth", "dispersion", "energy")



class TetrahedralSource(PhitsObject): # TODO: tetrahedral geometry
    name = "source"
    required = ["projectile", "energy"]
    positional = ["projectile", "energy"]
    optional = ["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
                "cell", "elevation", "azimuth", "dispersion"]
    ident_map = {"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
                 "weight": "wgt", "charge_override": "izst", "fissile": "ispfs",
                 "elevation": "dir", "azimuth": "phi", "dispersion": "dom", "energy": "e0"}
    value_map = {"neutrons": 2, True: 1}
    shape = ("s_type = 24", "projectile", "spin", "mask", "transform", "weight", "factor",
             "charge_override", "fissile", "mesh", "elevation", "azimuth", "dispersion", "energy")

class SurfaceSource(PhitsObject):
    name = "source"
    required = ["projectile", "energy", "surface", "cut"]
    positional = ["projectile", "energy", "surface", "cut"]
    optional = ["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile",
                "elevation", "azimuth", "dispersion"]
    ident_map = {"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl",
                 "weight": "wgt", "charge_override": "izst", "fissile": "ispfs",
                 "surface": "suf",
                 "elevation": "dir", "azimuth": "phi", "dispersion": "dom", "energy": "e0"}
    value_map = {"neutrons": 2, True: 1}
    shape = ("s_type = 26", "projectile", "spin", "mask", "transform", "weight", "factor",
             "charge_override", "fissile", "surface", "cut", "elevation", "azimuth", "dispersion", "energy")

# class Duct(PhitsObject):
#     name = "source"
#     required = ["wall", "dl0", "dl1", "dl2", "dpf", "drd"]
#     positional = ["wall", "dl0", "dl1", "dl2", "dpf", "drd"]
#     optional = ["dxw", "dyw"]
#     def __init__(self, *args, **kwargs):
