

class Rectangle(PhitsObject):
    name = "rectangle"
    syntax = {"width": ("tw", Real(), 0),
              "number": ("tn", Integer(), 1),
              "delta": ("td", Real(), 2),
              "center": ("t0", Real(), None)}
    grammar = r'start:  "t-type" "=" "1" "\n" @assign_among|self.inv_syntax| ~ 3..4'

    def assignment(self, tr):
        return (self.inv_syntax[tr[0]][0], self.inv_syntax[tr[0]][2].python(tr[1]))

    def start(self, tr):
        return Rectangle(**dict(tr))

class Gaussian(PhitsObject):
    name = "gaussian"
    syntax = {"fwhm": ("tw", Real(), 0),
              "number": ("tn", Integer(), 1),
              "delta": ("td", Real(), 2),
              "cutoff": ("tc", Real(), None)}
    grammar = r'start:  "t-type" "=" "2" "\n" @assign_among|self.inv_syntax| ~ 3..4'

    def assignment(self, tr):
        return (self.inv_syntax[tr[0]][0], self.inv_syntax[tr[0]][2].python(tr[1]))

    def start(self, tr):
        return Gaussian(**dict(tr))


class TimeDistributionBins(PhitsObject):
    name = "timedistbins"
    syntax = {"function": ("h(x)", Function(), 0),
              "n_bins": ("ll", Integer(), 1),
              "bounds": (("tg1", "tg2"), (Real(), Real()), 2),
              "particle_production": ("o-type", Array, None)}
    grammar = r'start: "t-type" "=" /3|4/ "\n" @assign_among|self.inv_syntax| ~ 4..5'


class TimeDistributionFunction(PhitsObject):
    name = "timedistfunction"
    syntax = {"function": ("h(x)", Function(), 0),
              "n_bins": ("ll", Integer(), 1),
              "bounds": (("tg1", "tg2"), (Real(), Real()), 2),
              "particle_production": ("o-type", Array, None)}
    grammar = r'start: "t-type" "=" /5|6/'

class TimeDistributionCustom(PhitsObject):
    name = "timedistcustom"
    syntax = {"function": ("h(x)", Function(), 0),
              "n_bins": ("ll", Integer(), 1),
              "bounds": (("tg1", "tg2"), (Real(), Real()), 2),
              "particle_production": ("o-type", Array, None)}
    grammar = r'start: "t-type" "=" /100/'

class TimeDistribution(PhitsObject):
    def __init__(self, function, n_bins, bounds, particle_production=None):
        

    def __init__(self, bins, weights, particle_production=None):
        assert len(bins) == len(weights), "The length of the bins and weights must be the same for a TimeDistribution."
        self.typ = 3
        self.bins = bins,
        self.weights = weights
        self.particle_production = particle_production



    mapping = {"function": ("h(x)", Function()),
               "n_bins": ("ll", None, Integer()),
               "adjust": }
    parser = r"""
    start: t0 | t3 | t4 | t5 | t6 | t100


    t3:  "t-type"  "="  "3"   "\n" assignment numbergrid
    t4:  "t-type"  "="  "4"  "\n" assignment numbergrid (assignment numberline?)?
    t5:  "t-type"  "="  "5"   "\n" assignment ~ 4
    t6:  "t-type"  "="  "6"  "\n" assignment ~ 4 (assignment numberline?)?
    t100:  "t-type"  "="  "100"  "\n" assignment ~ 2
    assignment: IDENTIFIER "=" computation "\n"
    """

    
    def t4(self, tr):
        if len(tr) <= 3:
            return TimeDistribution(tr[1][0], tr[1][1])
        else:
            return TimeDistribution(tr[1][0], tr[1][1], particle_production=tr[1][])

    def t3(self, tr):
        return TimeDistribution(tr[1][0], tr[1][1]) # I have literally no idea what the docs' "format" means

    def t0(self, tr):
        if tr[0] == "0":
            return
        elif tr[0] == "1":
            return Gaussian(tr[1], tr[2], tr[3], tr[4])
        elif tr[0] == "2":
            return Rectangle(tr[1], tr[2], tr[3], tr[4])



class AngleDistribution(PhitsObject):
    _parser = r"""
    atype: a1 | a4 | a5 | a6

    a1:  "a-type"  "="  ("1" | "11")   "\n" assignment numbergrid
    a4:  "a-type"  "="  ("4" | "14")   "\n" assignment numbergrid assignment numbergrid
    a5:  "a-type"  "="  ("5" | "15")   "\n" assignment ~ 4
    a6:  "a-type"  "="  /6|16/   "\n" assignment ~ 2 (assignment numbergrid?)?
    """

class EnergyDistribution(PhitsObject):
    _parser = r"""
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
    """
