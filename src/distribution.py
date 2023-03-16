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
    syntax = {"bins": (None, None, 0),
              "particle_production": (None, None, None)}

    shape = lambda self: ("t-type = 4" if self.particle_generation else "t-type = 3",
                          f"ntt = {len(self.bins)}",
                          "\n".join([" ".join(i) for i in self.bins]),
                          "o-type = 1\n" + " ".join(self.particle_production) if self.particle_production else "")


class AngleDistribution(PhitsObject):
    name = "source"
    syntax = {"bins": (None, None, 0),
              "unit": ("a-type", FinBij({"cos": 1, "degree": 11}), None),
              "particle_production": ("q-type", Choice10(), None)}

    shape = lambda self: ("a-type = 14" if self.unit == "degree" else "a-type = 4",
                          f"na = {len(self.bins)}",
                          "\n".join([" ".join(i) for i in self.bins]),
                          "q-type = 1\n" + " ".join(self.particle_production) if self.particle_production else "q-type = 0")


class EnergyDistribution(PhitsObject):
    name = "source"
    def __init__(self, bins, adjust="particles", units="MeV", normalize="1/MeV", particle_production=None):
        self.bins = bins
        self.adjust = adjust
        self.units = units
        self.normalize = normalize
        self.particle_production = particle_production

    shape = lambda self: (("e-type = 22" if self.units == "MeV" else "e-type = 32") if self.adjust == "particles" else \
                          ("e-type = 23" if self.units == "MeV" else "e-type = 33"),
                          f"ne = {-1 * len(self.bins)}" if self.normalize == "1/Lethargy" else f"ne = {len(self.bins)}",
                          "\n".join([" ".join(i) for i in self.bins]),
                          "p-type = 1\n" + " ".join(self.particle_production) if self.particle_production else "p-type = 0")
