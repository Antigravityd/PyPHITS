# Sources will be specified by an iterable whose elements are a tuple  (Source(), weight) where the last two are optional.
# totfact and iscorr may be specified in the general parameters section.

# Currently, no a-type mesh support for the elevation angles.
from base import *
from collections import namedtuple

common = {"type": ("s-type", None, finbij({})),
          "projectile": ("proj", None, finbij({})),
          "spin": (("sx", "sy", "sz"), (0, 0, 0), (posreal, posreal, posreal)),
          "mask": (("reg", "ntmax"), ("all", 1000), (IsA(Cell, index=True), posint)),
          "transform": ("trcl", "none", IsA(Transform, index=True)),
          "weight": ("wgt", 1.0, posreal),
          "charge_override": ("izst", "none", posreal),
          "counter_start": (("cnt(1)", "cnt(2)", "cnt(3)"), (0, 0, 0), (posint, posint, posint)),
          "fissile": ("ispfs", False, finbij({False: 0, "fissions": 1, "neutrons": 2}))
          # ibatch?
          }

semi_common = {"elevation": ("dir", 1.0, oneof(posreal, finbij({"isotropic": "all"}), isa(AngleDistribution, "data"))) #TODO:isa, oneof
               "azimuth": ("phi", "random", posreal), # TODO: random is sus
               "dispersion": ("dom", 0.0, oneof(posreal, finbij({"cos^2": -1})))
               "mono_energy": ("e0", None, posreal),
               "spectrum": ("e-type", None, isa(EnergyDistribution))}
# TODO: figure out how to deduce that exactly one of mono_energy, e-type is a required argument.



# TODO: implement parser macros and generation of transformers in PhitsObject methods
class Cylindrical(PhitsObject):
    mapping = common | {"center": (("x0", "y0"), (0.0, 0.0), (posreal, posreal)),
                        "zbounds": (("z0", "z1"), (0.0, 0.0), (posreal, posreal)),
                        "radius": ("r0", 0.0, posreal),
                        "cutout_radius": ("r1", 0.0, posreal)} | semi_common
    _parser = r"""
    start: "s-type" "=" "1" "\n" (assignment)+

    assignment: normal | etype | dir


    assignment: IDENTIFIER "=" computation "\n"

    normal: @assign_among([i for i in self.mapping.keys() if i not in ["spectrum", "elevation"]]) "\n"

    etype: @parseof(EnergyDistribution) // _LINEBRK?

    dir: "dir" "=" computation | "dir" "=" /all/ | @parseof(AngleDistribution) // _LINEBRK?
    """

    def dir(s, tr):
        if isinstance(tr[0], AngleDistribution):
            return

class Rectangular(PhitsObject):
    mapping = common | {"xbounds": (("x0", "x1"), (0.0, 0.0), (posreal, posreal)),
                        "ybounds": (("x0", "x1"), (0.0, 0.0), (posreal, posreal)),
                        "zbounds": (("x0", "x1"), (0.0, 0.0), (posreal, posreal))} | semi_common


class Source(PhitsObject): # currently no support for cnt(i) or ibatch common parameters

    mapping = {"type": ("s-type", None, finbij({})),
               "projectile": ("proj", None, finbij({})),
               "spin": (("sx", "sy", "sz"), (0, 0, 0), (posreal, posreal, posreal)),
               "mask": (("reg", "ntmax"), ("all", 1000), (isa(Cell), posint)),
               "transform": ("trcl", "none", isa(Transform)),
               "weight": ("wgt", 1.0, posreal),
               "charge_override": ("izst", "none", posreal),
               "counter_start": (("cnt(1)", "cnt(2)", "cnt(3)"), (0, 0, 0), (posint, posint, posint)),
               "fissile": ("ispfs", False, finbij({False: 0, "fissions": 1, "neutrons": 2}))
               # ibatch?
               "global_scaling":
               }
    parser = r"""
    %ignore SPACE
    start: assignment* (source | ("<source>" "=" computation source)+) assignment*

    source: "s-type" "=" POSINT (assignment | etype | ttype | atype)*


    etype: e1 | e2 | e3 | e4 | e5 | e6 | e7 | e20 | e25 | e28

    e1:  "e-type"  "="  ("1" | "8" | "11" | "18" | "21" | "22" | "31" | "32")  "\n" assignment numbergrid
    e4:  "e-type"  "="  ("4" | "9" | "14" | "19" | "23" | "24" | "33" | "34")   "\n" assignment numbergrid assignment numbergrid
    e2:  "e-type"  "="  ("2" | "12")   "\n" assignment ~ 4
    e3:  "e-type"  "="  ("3")   "\n" assignment ~
    e5:  "e-type"  "="  ("5" | "15")   "\n" assignment ~ 4
    e6:  "e-type"  "="  ("6" | "16")   "\n" assignment ~ 5 numbergrid
    e7:  "e-type"  "="  ("7")   "\n" assignment ~ 6 numbergrid
    e20:  "e-type"  "="  ("20")   "\n" assignment
    e25:  "e-type"  "="  ("25" | "26")   "\n" assignment ~ 14
    e28:  "e-type"  "="  ("28" | "29")   "\n" assignment ~ 7

    atype: a1 | a4 | a5 | a6

    a1:  "a-type"  "="  ("1" | "11")   "\n" assignment numbergrid
    a4:  "a-type"  "="  ("4" | "14")   "\n" assignment numbergrid assignment numbergrid
    a5:  "a-type"  "="  ("5" | "15")   "\n" assignment ~ 4
    a6:  "a-type"  "="  /6|16/   "\n" assignment ~ 2 (assignment numbergrid?)?

    ttype: t0 | t3 | t4 | t5 | t6  | t100

    t0:  "t-type"  "="  /0|1|2/   "\n" assignment ~ 5
    t3:  "t-type"  "="  "3"   "\n" assignment numbergrid
    t4:  "t-type"  "="  "4"  "\n" assignment numbergrid (assignment numbergrid?)?
    t5:  "t-type"  "="  "5"   "\n" assignment ~ 4
    t6:  "t-type"  "="  "6")  "\n" assignment ~ 4 (assignment numbergrid?)?
    t100:  "t-type"  "="  "100"  "\n" assignment ~ 2
    assignment: IDENTIFIER "=" computation "\n"
    """

    def assignment(self, ass):
        return (ass[0], ass[1])

    def t100(self, sec):
        raise NotImplementedError("Can't currently process things that require external files.")

    def t6(self, sec):
        weight_func = sec[0][1]
        n_weight_groups = sec[1][1]
        bounds = (sec[2][1], sec[3][1])
        generation_option = sec[4][1]
        if generation_option == 0:
            return ("time_distribution", {"function": weight_func, "n_bins": n_weight_groups, "bounds": bounds, "adjust": "weights"})
        else:
            return ("time_distribution", {"function": weight_func, "n_bins": n_weight_groups, "bounds": bounds, "adjust": "both",
                                          "particles_per_bin": sec[5]})

    def t5(self, sec):
        weight_func = sec[0][1]
        n_weight_groups = sec[1][1]
        bounds = (sec[2][1], sec[3][1])
        return ("time_distribution", {"function": weight_func, "n_bins": n_weight_groups, "bounds": bounds, "adjust": "particles"})

    def t4(self, sec):
        n_time_bins = sec[0][1]
        weights = sec[1]
        generation_option = sec[2][1]
        if generation_option == 0:
            return ("time_distribution", {"bins": weights, "adjust": "weights"})
        else:
            return ("time_distribution", {"bins": weights, "adjust": "weights", "particles_per_bin": sec[3]})

    def t3(self, sec):
        n_time_bins = sec[0][1]
        weights = sec[1]
        return ("time_distribution", {"bins": weights, "adjust": "particles"})

    def t0(self, sec):
        typ = sec[0]
        center = sec[1]
        width = sec[2]
        number = sec[3]
        delta = sec[4]
        cutoff = sec[5]
        if typ == 0:
            pass
        elif typ == 1:
            # TODO: update above to be tuples, and think about using namedtuple()
            return ("time_distribution", ("rectangle", [center, number, delta]))
        elif typ == 2:
            return ("time_distribution", ("gaussian", [center, number, delta, cutoff]))

    def a6(self, sec):
        unit = sec[0]
        func = sec[1][1]
        n_groups = sec[2][1]
        if len(sec) in [3, 4]:
            return ("angle_distribution", {"function": func, "n_groups": n_groups, "units": "degree" if unit == "6" else "degree",
                                           "adjust": "weights"})
        else:
            return ("angle_distribution", {"function": func, "n_groups": n_groups, "units": "degree" if unit == "6" else "degree",
                                           "adjust": "both", "particles_per_bin": sec[4]})

    def a5(self, sec):



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
