import sys
from base import *

_tester = Nuclide()
def _decomposition(composition):
    """Turns a composition list into PHITS input form."""
    r = ""
    for nuc, ratio in composition:
        conv = _tester.phits(nuc)
        if callable(conv):
            raise conv("composition")
        else:
            r += f"{conv} {ratio} "
    return r

class DataMax(PhitsObject):
    """Given a material, sets the maximum energy for an interaction between particles and nucleus in the material."""
    name = "data_max"
    syntax = {"particles": ("part", List(FinBij({"proton": "proton", "neutron": "neutron"}), unique=True), 0),
              "nucleus": (None, Nuclide(fake=True), 1),
              "max_energy": (None, PosReal(), 2)}
    superobjects = ["material"]
    prelude = ("particles", ("mat", "'nucleus", "dmax"))
    # Manual lies about accepting the usual syntax for nuclides
    shape = lambda self: (("material", "nucleus", "max_energy"),)
    group_by = lambda self: (self.particles,)
    separator = lambda self: self.section_title()
    max_groups = 6

class MatNameColor(PhitsObject):
    name = "mat_name_color"
    syntax = {"mat_name": (None, Text(), 0),
              "size": (None, PosReal(), 1),
              "color": (None, Color(), 2)} # TODO: color

    superobjects = ["material"]
    prelude = (("mat", "name", "\\size", "\\color"),)
    shape = (("material", "mat_name", "size", "color"),)


# TODO: temporarily setting the composition to be integer-only to cut down on length.
# TODO: also temporarily setting JENDL4Nuclide() to enable testing; it looks like Parameters() can change this error behavior,
# so should be set back to Nuclide() in production
# TODO: Nuclide-by-nuclide library setting
class Material(PhitsObject): # Composition is a list of pairs of (<element name string>, <ratio>) e.g. ("8Li", 0.5)
    name = "material"
    syntax = {"composition": (None, List(Tuple(JENDL4Nuclide(), PosInt()), unique_by=lambda x: kf_encode(x[0])), 0),
              "gas": ("GAS", Choice10(), None),
              "electron_step": ("ESTEP", PosInt(), None), # TODO: check integer right
              "neutron_lib": ("NLIB", LibraryID(), None),
              "photon_lib": ("PLIB", LibraryID(), None),
              "electron_lib": ("ELIB", LibraryID(), None),
              "proton_lib": ("HLIB", LibraryID(), None),
              "conductor": ("COND", FinBij({False: -1, True: 1}), None),
              "thermal_lib": (None, ThermalLib(), None),
              "chemical": ("chem", List(Tuple(Chemical(), PosInt())), None),
              "data_max": (None, IsA(DataMax, index=True), None),
              # "mat_time_change": (None, IsA(MatTimeChange, index=True), None),
              "mat_name_color": (None, IsA(MatNameColor, index=True), None)
              }
    subobjects = ["data_max", "mat_time_change", "mat_name_color"]
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



# TODO: irremovable circularity
# class MatTimeChange(PhitsObject):
#     name = "mat_time_change"
#     syntax = {"time": (None, PosReal(), 0, "non"),
#               "new": (None, IsA(Material, index=True), 1, "non"),
#               "old": (None, IsA(Material, index=True), None)}
#     prelude = (("mat", "'time", "change"))
#     shape = (("old", "time", "new"))



__pdoc__ = dict()
__pdoc__["builds"] = False
__pdoc__["slices"] = False
for name, cl in list(sys.modules[__name__].__dict__.items()):
    if type(cl) == type and issubclass(cl, PhitsObject) and cl != PhitsObject:
        __pdoc__[cl.__name__] = cl.__doc__ + cl.syntax_desc() if cl.__doc__ else cl.syntax_desc()
