from base import *

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

    shape = lambda self: ("t-type = 4" if self.particle_generation else "t-type = 3",
                          f"ntt = {len(self.bins)}",
                          "\n".join([" ".join(i) for i in self.bins]),
                          "o-type = 1\n" + " ".join(self.particle_production) if self.particle_production else "")


class AngleDistribution(PhitsObject):
    name = "source"
    syntax = {"bins": (None, List(PosReal()), 0),
              "unit": ("a-type", FinBij({"cos": 1, "degree": 11}), None),
              "particle_production": (None, List(PosReal()), None)}

    shape = lambda self: ("a-type = 14" if self.unit == "degree" else "a-type = 4",
                          f"na = {len(self.bins)}",
                          "\n".join([" ".join(i) for i in self.bins]),
                          "q-type = 1\n" + " ".join(self.particle_production) if self.particle_production else "q-type = 0")

class EnergyDistribution(PhitsObject):
    name = "source"
    syntax = {"bins": (None, List(PosReal()), None),
              "adjust": (None, FinBij({"particles": "particles", "weights": "weights"}), None),
              "units": (None, FinBij({"MeV": "MeV", "Angstrom": "Angstrom"}), None),
              "normalize": (None, FinBij({"1/Lethargy": -1, "1/MeV": 1}), None), # TODO: check 1/MeV
              "particle_production": (None, List(PosReal()), None)
              }

    shape = lambda self: (("e-type = 22" if self.units == "MeV" else "e-type = 32") if self.adjust == "particles" else \
                          ("e-type = 23" if self.units == "MeV" else "e-type = 33"),
                          f"ne = {self.normalize * len(self.bins)}",
                          "\n".join([" ".join(i) for i in self.bins]),
                          "p-type = 1\n" + " ".join(self.particle_production) if self.particle_production else "p-type = 0")

class GaussianEnergy(PhitsObject):
    name = "source"
    syntax = {"units": ("e-type", FinBij({"MeV": 2, "Angstrom": 12}), None, "MeV"),
              "center": ("eg0", Real(), None),
              "fwhm": ("eg1", PosReal(), None),
              "cutoffs": (("eg2", "eg3"), (PosReal(), PosReal()), None)}
    shape = lambda self: ("e-type = 2" if self.units == "MeV" else "e-type = 12", "center", "fwhm", "cutoffs")


class Maxwellian(PhitsObject):
    name = "source"
    syntax = {"particle_production": (None, List(PosReal()), None),
              "bins": (None, PosInt(), None),
              "interpolation": (None, FinBij({"linear": 1, "log": -1}), None),
              "temperature": ("et0", Real(), None),
              "cutoffs": (("et1", "et2"), (PosReal(), PosReal()), None),
              "power_index": ("et3", PosReal(), None)}

    shape = lambda self: ("e-type = 7" if self.particle_production else "e-type = 3",
                          f"nm = {str(self.interpolation * self.bins)}",
                          "termperature", "cutoffs", "power_index"
                          "p-type = 1\n" + " ".join(self.particle_production) if self.particle_production else "")

class Radioisotope(PhitsObject):
    name = "source"
    syntax = {"adjust": (None, FinBij({"particles": "particles", "weights": "weights"}), None),
              "isotopes": (None, List(Tuple(Nucleide(), PosReal())),  0),
              "decay": ("dtime", Real(), None), # TODO
              "min_activity": ("actlow", PosReal(), None),
              "normalize_per": ("norm", FinBij({"second": 0, "source": 1}), None),
              "electron_production": ("iaugers", FinBij({"all": 0, "beta": 1, "auger/IC": 2})),
              "annihilation_photons": ("iannih", Choice10(c_style=True), None),
              "photon_production": ("icharctx", FinBij({"all": 0, "gamma": 1, "xray": 2}), None)}

    shape = lambda self: ("e-type = 28" if self.adjust == "particles" else "e-type = 29",
                          f"ni = {len(self.isotopes)}",
                          "isotopes", "decay", "min_activity", "normalize_per", "electron_production", "annihilation_photons",
                          "photon_production")

# TODO: cosmic ray
