

"""Anything that looks like `a-type = ...` in a PHITS .inp file.
Generally represents some sort of distribution.
"""

from base import *
import sys




class Rectangle(PhitsObject):
    name = "source"
    syntax = {"width": ("tw", Real(), 0),
              "number": ("tn", PosInt(), 1),
              "delta": ("td", Real(), 2),
              "center": ("t0", Real(), None)}

    shape = ("t-type = 1", "center", "width", "number", "delta")


class Gaussian(PhitsObject):
    name = "source"
    syntax = {"fwhm": ("tw", Real(), 0),
              "number": ("tn", Integer(), 1),
              "delta": ("td", Real(), 2),
              "cutoff": ("tc", Real(), None)}

    shape = ("t-type = 2", "fwhm", "number", "delta", "cutoff")




class TimeDistribution(PhitsObject):
    name = "source"
    syntax = {"bins": (None, List(PosReal()), 0),
              "particle_production": (None, List(PosReal()), None)}

    shape = lambda slf: ("t-type = 4" if slf.particle_generation else "t-type = 3",
                         f"ntt = {len(slf.bins)}",
                         "\n".join(" ".join(str(j) for j in i) for i in slf.bins),
                         "o-type = 1\n" + " ".join((str(i) for i in slf.particle_production)) if slf.particle_production else "")


class AngleDistribution(PhitsObject):
    name = "source"
    syntax = {"bins": (None, List(Tuple(PosReal(), PosReal())), 0),
              "last_bin": (None, PosReal(), 1),
              "unit": ("a-type", FinBij({"cos": 1, "degree": 11}), None),
              "particle_production": (None, List(PosReal()), None)}

    shape = lambda self: ("a-type = 14" if self.unit == "degree" else "a-type = 4",
                          f"na = {len(self.bins)}",
                          "\n".join(" ".join(str(i) for i in j) for j in self.bins) + f"\n{self.last_bin}",
                          "q-type = 1\n" + " ".join((str(i) for i in self.particle_production)) if self.particle_production \
                          else "q-type = 0")

    def restrictions(self):
        if len(self.bins) != len(self.particle_production):
            raise ValueError("For EnergyDistribution: len(bins) must equal len(particle_production).")



class EnergyDistribution(PhitsObject):
    name = "source"
    syntax = {"bins": (None, List(Tuple(PosReal(), PosReal(), PosReal())), None),
              "adjust": (None, FinBij({"particles": "particles", "weights": "weights"}), None),
              "units": (None, FinBij({"MeV": "MeV", "Angstrom": "Angstrom"}), None),
              "normalize": (None, FinBij({"1/Lethargy": -1, "1/MeV": 1}), None), # TODO: check 1/MeV
              "particle_production": (None, List(PosReal()), None)
              }

    shape = lambda slf: (("e-type = 22" if slf.units == "MeV" else "e-type = 32") if slf.adjust == "particles" else \
                          ("e-type = 23" if slf.units == "MeV" else "e-type = 33"),
                          f"ne = -{len(slf.bins)}" if slf.normalize == "1/Lethargy" else f"ne = {len(slf.bins)}",
                          "\n".join(" ".join(str(j) for j in i) for i in slf.bins),
                          "p-type = 1\n" + " ".join((str(i) for i in slf.particle_production)) if slf.particle_production else "p-type = 0")
    def restrictions(self):
        if len(self.bins) != len(self.particle_production):
            raise ValueError("For EnergyDistribution: len(bins) must equal len(particle_production).")

class GaussianEnergy(PhitsObject):
    name = "source"
    syntax = {"units": ("e-type", FinBij({"MeV": 2, "Angstrom": 12}), None, "MeV"),
              "center": ("eg0", Real(), None),
              "fwhm": ("eg1", PosReal(), None),
              "cutoffs": (("eg2", "eg3"), (PosReal(), PosReal()), None)}
    shape = lambda slf: ("e-type = 2" if slf.units == "MeV" else "e-type = 12", "center", "fwhm", "cutoffs")


class Maxwellian(PhitsObject):
    name = "source"
    syntax = {"particle_production": (None, List(PosReal()), None),
              "bins": (None, PosInt(), None),
              "interpolation": (None, FinBij({"linear": 1, "log": -1}), None),
              "temperature": ("et0", Real(), None),
              "cutoffs": (("et1", "et2"), (PosReal(), PosReal()), None),
              "power_index": ("et3", PosReal(), None)}

    shape = lambda slf: ("e-type = 7" if slf.particle_production else "e-type = 3",
                          f"nm = {str(slf.interpolation * slf.bins)}",
                          "termperature", "cutoffs", "power_index"
                          "p-type = 1\n" + " ".join((str(i) for i in slf.particle_production)) if slf.particle_production else "")

class Radioisotope(PhitsObject):
    name = "source"
    syntax = {"adjust": (None, FinBij({"particles": "particles", "weights": "weights"}), None),
              "isotopes": (None, List(Tuple(Nuclide(), PosReal())),  0),
              "decay": ("dtime", Real(), None), # TODO
              "min_activity": ("actlow", PosReal(), None),
              "normalize_per": ("norm", FinBij({"second": 0, "source": 1}), None),
              "electron_production": ("iaugers", FinBij({"all": 0, "beta": 1, "auger/IC": 2}), None),
              "annihilation_photons": ("iannih", Choice10(c_style=True), None),
              "photon_production": ("icharctx", FinBij({"all": 0, "gamma": 1, "xray": 2}), None)}

    shape = lambda slf: ("e-type = 28" if slf.adjust == "particles" else "e-type = 29",
                          f"ni = {len(slf.isotopes)}",
                          "isotopes", "decay", "min_activity", "normalize_per", "electron_production", "annihilation_photons",
                          "photon_production")

# TODO: cosmic ray
__pdoc__ = dict()
__pdoc__["builds"] = False
__pdoc__["slices"] = False
for name, cl in list(sys.modules[__name__].__dict__.items()):
    if type(cl) == type and issubclass(cl, PhitsObject) and cl != PhitsObject: # and cl.__module__ == __name__
        __pdoc__[cl.__name__] = cl.__doc__ + cl.syntax_desc() if cl.__doc__ else cl.syntax_desc()
