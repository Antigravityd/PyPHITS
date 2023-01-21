2
from base import *


class MatTimeChange(PhitsObject):
    name = "mat_time_change"
    required = ["time", "new"]
    positional = ["time", "new"]
    optional = ["old"]
    shape = (("old", "time", "new"))
    prelude = (("mat", "\\time", "change"))
    nones = {"time": "non", "change": "non"}



class DataMax(PhitsObject): # requires special handling in make_input
    name = "data_max"
    required = ["particles", "nucleus", "threshold"]
    positional = ["particles", "nucleus", "threshold"]
    optional = ["material"]
    shape = (("material", "nucleus", "threshold"))
    prelude = ("particles", (("mat", "\\nucleus", "dmax")))
    group_by = lambda self: self.particles
    max_groups = 6
    separator = lambda self: self.section_title()
    ident_map = {"particles": "part"}


class MatNameColor(PhitsObject):
    name = "mat_name_color"
    required = ["name", "size", "color"]
    positional =  ["name", "size", "color"]
    optional = ["material"]
    shape = (("material", "name", "size", "color"))
    prelude = (("mat", "\\name", "\\size", "\\color"))



# TODO: how the libraries work isn't well-documented. Is there a single library set for the whole material, or
# does one set a library after each element of the compositon? Can the thermal neutron library be set anywhere?
class Material(PhitsObject): # Composition is a list of pairs of (<element name string>, <ratio>) e.g. ("8Li", 0.5)
    name = "material"
    required = ["composition"]
    positional = ["composition"]
    optional = ["chemical", "time_change", "data_max", "mat_name_color", "condensed", "conductive", "electron_step",
                "neutron_lib", "proton_lib", "electron_lib", "photon_lib", "thermal_lib"]
    shape = (lambda self: f"MAT[{self.index}]",
             (lambda self: "".join(map(lambda tup: f"{tup[0]} {tup[1]} ", self.composition))), "chemical", "condensed",
             "conductive", "electron_step", "neutron_lib", "proton_lib","electron_lib", "photon_lib",
             lambda self: f"MT{self.index} {self.thermal_lib}" if self.thermal_lib is not None else "")
    subobjects = ["time_change", "data_max", "mat_name_color"]
    ident_map = {"condensed": "GAS", "conductive": "COND", "electron_step": "ESTEP", "neutron_lib": "NLIB",
                 "proton_lib": "HLIB", "electron_lib": "HLIB", "photon_lib": "PLIB", "chemical": "chem"}
