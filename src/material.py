import sys
from base import *

_tester = Nuclide()
def _decomposition(composition):
    r = ""
    for nuc, ratio in composition:
        conv = _tester.phits(nuc)
        if callable(conv):
            raise conv("composition")
        else:
            r += f"{conv} {ratio} "
    return r



# TODO: how the libraries work isn't well-documented. Is there a single library set for the whole material, or
# does one set a library after each element of the compositon? Can the thermal neutron library be set anywhere?
# TODO: temporarily setting the composition to be integer-only to cut down on length.
class Material(PhitsObject): # Composition is a list of pairs of (<element name string>, <ratio>) e.g. ("8Li", 0.5)
    name = "material"
    syntax = {"composition": (None, List(Tuple(JENDL4Nuclide(), PosInt())), 0),
              "gas": ("GAS", Choice10(), None),
              "electron_step": ("ESTEP", PosInt(), None), # TODO: check integer right
              "neutron_lib": ("NLIB", LibraryID(), None),
              "photon_lib": ("PLIB", LibraryID(), None),
              "electron_lib": ("ELIB", LibraryID(), None),
              "proton_lib": ("HLIB", LibraryID(), None),
              "conductor": ("COND", FinBij({False: -1, True: 1}), None),
              "thermal_lib": (None, ThermalLib(), None),
              "chemical": ("chem", List(Tuple(Chemical(), PosInt())), None),
              # "mat_time_change": (None, IsA(MatTimeChange, index=True), None),
              }
    shape = lambda self: (f"MAT[{self.index}]",
                          _decomposition(self.composition),
                          "gas", "electron_step", # "neutron_lib", "photon_lib", "electron_lib", "proton_lib",
                          "conductor",
                          "chem = " + " ".join(ch + " " + str(den) for ch, den in self.chemical) \
                          if self.chemical is not None else "",
                          f"MT{self.index} " + self.syntax["thermal_lib"][1].phits(self.thermal_lib) \
                          if self.thermal_lib is not None else ""
                          )

    def restrictions(self):
        if any(map(lambda x: int(str(_tester.phits(x[0]))[:-3]) > 92, self.composition)) and \
           (not hasattr(self, "parameters") or not hasattr(self.parameters, "stoping_model") or \
            "ATIMA" in self.parameters.stopping_model):
            raise ValueError("Material cannot have nuclei with Z > 92 and ATIMA set at the same time;"
                             " please pass stopping_model=SPAR+NTMC to the material.")

        if len(set(map(lambda x: x[0], self.composition))) < len(self.composition):
            raise ValueError(f"Material cannot have duplicate nuclei; got {self.composition}")




# TODO: irremovable circularity
# class MatTimeChange(PhitsObject):
#     name = "mat_time_change"
#     syntax = {"time": (None, PosReal(), 0, "non"),
#               "new": (None, IsA(Material, index=True), 1, "non"),
#               "old": (None, IsA(Material, index=True), None)}
#     prelude = (("mat", "'time", "change"))
#     shape = (("old", "time", "new"))

class MatNameColor(PhitsObject):
    name = "mat_name_color"
    syntax = {"name": (None, Text(), 0),
              "size": (None, PosReal(), 1),
              "color": (None, Color(), 2)} # TODO: color

    superobjects = ["material"]
    prelude = (("mat", "\\name", "\\size", "\\color"),)
    shape = (("material", "name", "size", "color"),)


__pdoc__ = dict()
__pdoc__["builds"] = False
__pdoc__["slices"] = False
for name, cl in list(sys.modules[__name__].__dict__.items()):
    if type(cl) == type and issubclass(cl, PhitsObject) and cl != PhitsObject:
        __pdoc__[cl.__name__] = cl.__doc__ + cl.syntax_desc() if cl.__doc__ else cl.syntax_desc()
