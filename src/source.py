# Sources will be specified by an iterable whose elements are a tuple  (Source(), weight) where the last two are optional.
# totfact and iscorr may be specified in the general parameters section.

# Currently, no a-type mesh support for the elevation angles.
from base import *

class Source(PhitsObject): # currently no support for cnt(i) or ibatch common parameters
    parser = r"""
    start: assignment* (source | ("<source>" "=" computation source)+) assignment*

    source: "s-type" "=" POSINT (assignment | etype | ttype | atype)*

    assignment: IDENTIFIER "=" computation


    etype: e1 | e2 | e3 | e4 | e5 | e6 | e7 | e20 | e25 | e28

    e1: " "* "e-type"i " "* "=" " "* /1|8|11|18|21|22|31|32/ " "* "\n" assignment numbergrid
    e4: " "* "e-type"i " "* "=" " "* /4|9|14|19|23|24|33|34/  " "* "\n" assignment numbergrid assignment numbergrid
    e2: " "* "e-type"i " "* "=" " "* /2|12/  " "* "\n" assignment ~ 4
    e3: " "* "e-type"i " "* "=" " "* /3/  " "* "\n" assignment ~
    e5: " "* "e-type"i " "* "=" " "* /5|15/  " "* "\n" assignment ~ 4
    e6: " "* "e-type"i " "* "=" " "* /6|16/  " "* "\n" assignment ~ 5 numbergrid
    e7: " "* "e-type"i " "* "=" " "* /7/  " "* "\n" assignment ~ 6 numbergrid
    e20: " "* "e-type"i " "* "=" " "* /20/  " "* "\n" assignment
    e25: " "* "e-type"i " "* "=" " "* /25|26/  " "* "\n" assignment ~ 14
    e28: " "* "e-type"i " "* "=" " "* /28|29/  " "* "\n" assignment ~ 7

    atype: p1 | p4 | p5 | p6

    p1: " "* "a-type"i " "* "=" " "* /1|11/  " "* "\n" assignment numbergrid
    p4: " "* "a-type"i " "* "=" " "* /4|14/  " "* "\n" assignment numbergrid assignment numbergrid
    p5: " "* "a-type"i " "* "=" " "* /5|15/  " "* "\n" assignment ~ 4
    p6: " "* "a-type"i " "* "=" " "* /6|16/  " "* "\n" assignment ~ 3 numbergrid

    ttype: t0 | t3 | t4 | t5 | t6  | t100

    t0: " "* "t-type"i " "* "=" " "* /0|1|2/  " "* "\n" assignment ~ 5
    t3: " "* "t-type"i " "* "=" " "* /3/  " "* "\n" assignment numbergrid
    t4: " "* "t-type"i " "* "=" " "* /0|1|2/  " "* "\n" assignment numbergrid assignment numbergrid
    t5: " "* "t-type"i " "* "=" " "* /5/  " "* "\n" assignment ~ 4
    t6: " "* "t-type"i " "* "=" " "* /6/  " "* "\n" assignment ~ 5 numbergrid
    t100: " "* "t-type"i " "* "=" " "* /100/  " "* "\n" assignment ~ 2
    """

    def start(self, tree):

    def source

    def __init__(self, s_type, projectile, *, spin=(0, 0, 0), mask=(None, 1000),
                 transform=None, weight=1.0, factor=1.0, charge_override=None, fissile=False, **kwargs):
        super().__init__("source", **kwargs)
        self.s_type = s_type
        self.proj = projectile
        self.sx = spin[0]
        self.sy = spin[1]
        self.sz = spin[2]
        self.reg = mask[0]
        self.ntmax = mask[1]
        self.trcl = transform
        self.wgt = weight
        self.factor = factor
        self.izst = charge_override
        self.ispfs = 0 if not fissile else (2 if fissile == "neutrons" else 1)

    def definition(self):
        inp = ""
        for var, val in [(k, v) for k, v in self.__dict__.items() if k not in {"index", "name", "parameters"}]:
            if val is not None:
                if isinstance(val, PhitsObject):
                    inp += f"{var} = {val.index}\n"
                else:
                    var2 = var.replace("_", "-")
                    inp += f"{var2} = {val}\n"

        return inp

        

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
