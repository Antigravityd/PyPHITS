
from base import *

# TODO: how the libraries work isn't well-documented. Is there a single library set for the whole material, or
# does one set a library after each element of the compositon? Can the thermal neutron library be set anywhere?
class Material(PhitsObject): # Composition is a list of pairs of (<element name string>, <ratio>) e.g. ("8Li", 0.5)
    name = "material"
    syntax = {"composition": (None, List(Tuple(Element(), PosReal())), 0), # TODO: tuple, element
              "gas": ("GAS", Choice10(), None),
              "electron_step": ("ESTEP", Integer(), None), # TODO: check integer right
              "neutron_lib": ("NLIB", Integer(), None),
              "photon_lib": ("PLIB", Integer(), None),
              "electron_lib": ("ELIB", Integer(), None),
              "proton_lib": ("HLIB", Integer(), None),
              "conductor": ("COND", FinBij({False: -1, True: 1}), None),
              "thermal_lib": (None, Integer(), None),
              "chemical": ("chem", List(Tuple(Chemical(), PosReal())), None), # TODO: chemical
              "time_change": (None, IsA(MatTimeChange), None),
              "data_max": (None, IsA(DataMax), None)
              }
    shape = lambda self: (f"MAT[{self.index}]",
                          "".join(map(lambda tup: f"{tup[0]} {tup[1]} ", self.composition)),
                          "gas", "electron_step", "neutron_lib", "photon_lib", "electron_lib", "proton_lib", "conductor",
                          "thermal_lib", "chemical",
                          f"MT{self.index} {self.thermal_lib}" if self.thermal_lib is not None else "")

    subobjects = ["time_change", "data_max", "mat_name_color"]



class MatTimeChange(PhitsObject):
    name = "mat_time_change"
    syntax = {"time": (None, PosReal(), 0, "non"),
              "new": (None, IsA(Material, index=True), 1, "non"),
              "old": (None, IsA(Material, index=True), None)}
    prelude = (("mat", "'time", "change"))
    shape = (("old", "time", "new"))



class DataMax(PhitsObject): # requires special handling in make_input
    name = "data_max"
    syntax = {"particles": ("part", List(FinBij({"neutron": "neutron", "proton": "proton", "all": "all"})), 0), # TODO: particle
              "nucleus": (None, OneOf(Element(), FinBij({"all": "all"})), 1),
              "threshold": (None, PosReal(), 2),
              "material": (None, OneOf(IsA(Material), FinBij({"all": "all"})), None)
              }
    prelude = ("particles", (("mat", "'nucleus", "dmax")))
    shape = (("material", "nucleus", "threshold"))
    group_by = lambda self: self.particles
    max_groups = 6
    separator = lambda self: self.section_title()


# TODO: necessary?
# class MatNameColor(PhitsObject):
#     name = "mat_name_color"
#     required = ["name", "size", "color"]
#     positional =  ["name", "size", "color"]
#     optional = ["material"]
#     shape = (("material", "name", "size", "color"))
#     prelude = (("mat", "\\name", "\\size", "\\color"))



