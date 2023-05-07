from functools import partial
from hypothesis.strategies import *
import re
from scipy.stats import ortho_group
import numpy as np
from numpy.linalg import det

# TODO: think about if the python() methods are necessary



class ValSpec():
    def __init__(self, strat):
        self.strat = strat
    def __or__(self, other):
        # TODO: think about if copying is necessary
        if isinstance(self, OneOf) and isinstance(other, OneOf):
            self.choices.append(other.choices)
            return self
        elif isinstance(self, OneOf) and not isinstance(other, OneOf):
            self.choices.append(other)
            return self
        elif not isinstance(self, OneOf) and isinstance(other, OneOf):
            other.choices.append(self)
            return other
        elif not isinstance(self, OneOf) and not isinstance(other, OneOf):
            return OneOf(self, other)


class Choice10(ValSpec):
    def __init__(self, c_style=False, true=True, false=False):
        super().__init__(one_of(just(false), just(true)))
        self.c_style = c_style
        self.true = true
        self.false = false

    def phits(self, val):
        if val == self.true:
            return 0 if self.c_style else 1
        elif val == self.false:
            return 1 if self.c_style else 0
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be either True or False; got {va}."), val)

    def python(self, val):
        if val == 0:
                return self.true if self.c_style else self.false
        elif val == 1:
            return self.false if self.c_style else self.true
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be either 0 or 1; got {va}."), val)

    def description(self):
        return f"either {self.true} or {self.false}"


class Integer(ValSpec):
    def __init__(self):
        super().__init__(integers())
    def phits(self, val):
        if isinstance(val, int):
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be an integer; got {va}."), val)

    def python(self, val):
        if val % 1 == 0:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be an integer; got {va}."), val)
    def description(self):
        return "int"

class Real(ValSpec):
    def __init__(self):
        super().__init__(floats(allow_nan=False, allow_infinity=False, allow_subnormal=False))
    def phits(self, val):
        if isinstance(val, float) or isinstance(val, int):
            return float(val)
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a float; got {va}."), val)

    def python(self, val):
        return val

    def description(self):
        return "float"

class PosInt(ValSpec):
    def __init__(self):
        super().__init__(integers(min_value=1))

    def phits(self, val):
        if isinstance(val, int) and val > 0:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a positive integer; got {va}."), val)

    def python(self, val):
        if val > 0:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be positive; got {va}."), val)

    def description(self):
        return "int > zero"


class PosReal(ValSpec):
    def __init__(self):
        super().__init__(floats(min_value=0, exclude_min=True, allow_nan=False, allow_infinity=False, allow_subnormal=False))
    def phits(self, val):
        if (isinstance(val, float) or isinstance(val, int)) and val > 0:
            return float(val)
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a positive float; got {va}."), val)

    def python(self, val):
        if val > 0:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be positive; got {va}."), val)

    def description(self):
        return "float > 0"




class NegDisable(ValSpec):
    def __init__(self):
        super().__init__(one_of(none(), integers(min_value=0), floats(min_value=0, allow_nan=False, allow_infinity=False,
                                                                      allow_subnormal=False)))
    def phits(self, val):
        if isinstance(val, float) or isinstance(val, int) and val >= 0:
            return val
        elif val is None:
            return -1.0
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a positive integer or None; got {va}."), val)

    def python(self, val):
        if val > 0:
            return val
        else:
            return None

    def description(self):
        return "either None or a number"

class RealBetween(ValSpec):
    def __init__(self, start, stop):
        super().__init__(floats(min_value=start, max_value=stop, allow_nan=False, allow_infinity=False, allow_subnormal=False))
        self.start = start
        self.stop = stop

    def phits(self, val):
        if isinstance(val, float) and (self.start is None or val >= self.start) and (self.stop is None or val <= self.stop):
            return val
        else:
            return partial(lambda va, start, stop, var: ValueError(f"`{var}` must be a float between {start} and {stop}; got {va}."),
                           val, self.start, self.stop)


    def python(self, val):
        if val >= self.start and val <= self.stop:
            return val
        else:
            return partial(lambda va, start, stop, var: ValueError(f"`{var}` must be a float between {start} and {stop}; got {va}."),
                           val, self.start, self.stop)

    def description(self):
        return f"float between {self.start} and {self.stop}"

class Between(ValSpec):
    def __init__(self, start, stop):
        super().__init__(integers(min_value=start, max_value=stop))
        self.start = start
        self.stop = stop

    def phits(self, val):
        if isinstance(val, int) and val >= self.start and val <= self.stop:
            return val
        else:
            return partial(lambda va, start, stop, var: ValueError(f"`{var}` must be an integer between {start} and {stop}; got {va}."),
                           val, self.start, self.stop)


    def python(self, val):
        if val >= self.start and val <= self.stop:
            return val
        else:
            return partial(lambda va, start, stop, var: ValueError(f"`{var}` must be an integer between {start} and {stop}; got {va}."),
                           val, self.start, self.stop)

    def description(self):
        return f"integer between {self.start} and {self.stop}"



class ZeroSpecial(ValSpec):
    def __init__(self, zero):
        super().__init__(one_of(just(zero), integers()))
        self.zero = zero

    def phits(self, val):
        if isinstance(val, int):
            if val == self.zero:
                return 0
            else:
                return val
        else:
            return partial(lambda va, zero, var: ValueError(f"`{var}` must be an integer or {zero}; got {va}."), val, self.zero)

    def python(self, val):
        if val == 0:
            return self.zero
        else:
            return val

    def description(self):
        return f"either {self.zero} or an integer"

class FinBij(ValSpec):
    def __init__(self, dic):
        super().__init__(sampled_from(list(dic.keys())))
        self.dic = dic

    def phits(self, val):
        if val in self.dic:
            return self.dic[val]
        else:
            return partial(lambda va, keys, var: ValueError(f"`{var}` must be one of {keys}; got {va}."), val, list(self.dic.keys()))


    def python(self, val):
        rev = {v: k for k, v in self.dic.items()}
        if val in rev:
            return rev[val]
        else:
            return partial(lambda va, keys, var: ValueError(f"`{var}` must be one of {keys}; got {va}."), val, list(rev.keys()))

    def description(self):
        return f"one of the keys in {self.dic}, with the value being the corresponding PHITS value"


class IsA(ValSpec):
    def __init__(self, cls, index=False):
        req = list(map(lambda tup: tuples(*[i.strat for i in tup[1]]) if isinstance(tup[1], tuple) else tup[1].strat,
                       sorted([v for k, v in cls.syntax.items() if v[2] is not None], key=lambda tup: tup[2])))
        opt = dict(map(lambda tup: (tup[0], tuples(*[i.strat for i in tup[1][1]]) if isinstance(tup[1][1], tuple) else tup[1][1].strat),
                       [(k, v) for k, v in cls.syntax.items() if v[2] is None]))

        super().__init__(builds(cls, *req, **opt)) # TODO: sus
        self.cls = cls
        self.index = index

    def phits(self, val):
        if not isinstance(val, self.cls):
            return partial(lambda va, cls, var: ValueError(f"`{var}` must be an instance of {cls}; got {val}."), val, self.cls)

        if self.index:
            if self.index is None:
                breakpoint()
            return val.index
        else:
            return val.definition()

    def python(self, val):
        return val

    def description(self):
        return f"an instance of {self.cls.__name__}"



class List(ValSpec):
    def __init__(self, entr):
        super().__init__(lists(entr.strat, min_size=1))
        self.entr = entr

    def phits(self, val):
        return " ".join(map(lambda v: str(self.entr.phits(v)), val))

    def python(self, val):
        return " ".join(map(lambda v: str(self.entr.python(v)), val))

    def description(self):
        return f"a list of values, each of which is {self.entr.description()}"

class Tuple(ValSpec):
    def __init__(self, *entr):
        super().__init__(tuples(*[i.strat for i in entr]))
        self.entr = entr

    def phits(self, val):
        return tuple(map(lambda t: self.entr[t[0]].phits(t[1]), enumerate(val)))

    def python(self, val):
        return tuple(map(lambda t: self.entr[t[0]].python(t[1]), enumerate(val)))

    def description(self):
        j = ", "
        return f"({j.join((i.description() for i in self.entr))})"


class OneOf(ValSpec):
    def __init__(self, *args):
        assert all(map(lambda x: isinstance(x, ValSpec), args)), "All arguments to OneOf must be value specifications."
        super().__init__(one_of(*[i.strat for i in args]))
        self.choices = args

    def phits(self, val):
        return self.that_which_applies(val, "phits").phits(val)

    def python(self, val):
        return self.that_which_applies(val, "python").python(val)

    def that_which_applies(self, val, wh):
        def _applies(s, val, wh):
            if wh == "phits":
                r = s.phits(val)
                if callable(r) or isinstance(r, Exception):
                    return False
                else:
                    return True
            else:
                r = s.python(val)
                if callable(r) or isinstance(r, Exception):
                    return False
                else:
                    return True

        applicable = list(filter(lambda x: _applies(x, val, wh), self.choices))
        if len(applicable) == 0:
            breakpoint()
            raise ValueError("Empty OneOf value specification.")
        else:
            if len(applicable) != 1:
                breakpoint()
            assert len(applicable) == 1, "Ambiguous OneOf value specification."
            return applicable[0]

    def description(self):
        return "either " + ", ".join(map(lambda x: x.description(), self.choices[:-1])) + ", or " + self.choices[-1].description()

elements = {1: ("H", "Hydrogen"), 2: ("He", "Helium"), 3: ("Li", "Lithium"), 4: ("Be", "Beryllium"), 5: ("B", "Boron"),
            6: ("C", "Carbon"), 7: ("N", "Nitrogen"), 8: ("O", "Oxygen"), 9: ("F", "Fluorine"), 10: ("Ne", "Neon"),
            11: ("Na", "Sodium"), 12: ("Mg", "Magnesium"), 13: ("Al", "Aluminum"), 14: ("Si", "Silicon"),
            15: ("P", "Phosphorus"), 16: ("S", "Sulfur"), 17: ("Cl", "Chlorine"), 18: ("Ar", "Argon"), 19: ("K", "Potassium"),
            20: ("Ca", "Calcium"), 21: ("Sc", "Scandium"), 22: ("Ti", "Titanium"), 23: ("V", "Vanadium"),
            24: ("Cr", "Chromium"), 25: ("Mn", "Manganese"), 26: ("Fe", "Iron"), 27: ("Co", "Cobalt"), 28: ("Ni", "Nickel"),
            29: ("Cu", "Copper"), 30: ("Zn", "Zinc"), 31: ("Ga", "Gallium"), 32: ("Ge", "Germanium"), 33: ("As", "Arsenic"),
            34: ("Se", "Selenium"), 35: ("Br", "Bromine"), 36: ("Kr", "Krypton"), 37: ("Rb", "Rubidium"),
            38: ("Sr", "Strontium"), 39: ("Y", "Yttrium"), 40: ("Zr", "Zirconium"), 41: ("Nb", "Niobium"),
            42: ("Mo", "Molybdenum"), 43: ("Tc", "Technetium"), 44: ("Ru", "Ruthenium"), 45: ("Rh", "Rhodium"),
            46: ("Pd", "Palladium"), 47: ("Ag", "Silver"), 48: ("Cd", "Cadmium"), 49: ("In", "Indium"), 50: ("Sn", "Tin"),
            51: ("Sb", "Antimony"), 52: ("Te", "Tellurium"), 53: ("I", "Iodine"), 54: ("Xe", "Xenon"), 55: ("Cs", "Cesium"),
            56: ("Ba", "Barium"), 57: ("La", "Lanthanum"), 58: ("Ce", "Cerium"), 59: ("Pr", "Praseodymium"),
            60: ("Nd", "Neodymium"), 61: ("Pm", "Promethium"), 62: ("Sm", "Samarium"), 63: ("Eu", "Europium"),
            64: ("Gd", "Gadolinium"), 65: ("Tb", "Terbium"), 66: ("Dy", "Dysprosium"), 67: ("Ho", "Holmium"),
            68: ("Er", "Erbium"), 69: ("Tm", "Thulium"), 70: ("Yb", "Ytterbium"), 71: ("Lu", "Lutetium"), 72: ("Hf", "Hafnium"),
            73: ("Ta", "Tantalum"), 74: ("W", "Tungsten"), 75: ("Re", "Rhenium"), 76: ("Os", "Osmium"), 77: ("Ir", "Iridium"),
            78: ("Pt", "Platinum"), 79: ("Au", "Gold"), 80: ("Hg", "Mercury"), 81: ("Tl", "Thallium"), 82: ("Pb", "Lead"),
            83: ("Bi", "Bismuth"), 84: ("Po", "Polonium"), 85: ("At", "Astatine"), 86: ("Rn", "Radon"), 87: ("Fr", "Francium"),
            88: ("Ra", "Radium"), 89: ("Ac", "Actinium"), 90: ("Th", "Thorium"), 91: ("Pa", "Protactinium"),
            92: ("U", "Uranium"), 93: ("Np", "Neptunium"), 94: ("Pu", "Plutonium"), 95: ("Am", "Americium"),
            96: ("Cm", "Curium"), 97: ("Bk", "Berkelium"), 98: ("Cf", "Californium"), 99: ("Es", "Einsteinium"),
            100: ("Fm", "Fermium"), 101: ("Md", "Mendelevium"), 102: ("No", "Nobelium"), 103: ("Lr", "Lawrencium"),
            104: ("Rf", "Rutherfordium"), 105: ("Db", "Dubnium"), 106: ("Sg", "Seaborgium"), 107: ("Bh", "Bohrium"),
            108: ("Hs", "Hassium"), 109: ("Mt", "Meitnerium"), 110: ("Ds", "Darmstadtium"), 111: ("Rg", "Roentgenium"),
            112: ("Cn", "Copernicium"), 113: ("Nh", "Nihonium"), 114: ("Fl", "Flerovium"), 115: ("Mc", "Moscovium"),
            116: ("Lv", "Livermorium"), 117: ("Ts", "Tennessine"), 118: ("Og", "Oganesson")}

particles = {2212: "proton", 2112: "neutron", 211: "pion+", 111: "pion0", -211: "pion-", -13: "muon+", 13: "muon-",
             321: "kaon+", 311: "kaon0", -321: "kaon-", 11: "electron", -11: "positron", 22: "photon",
             12: 'e_neutrino', -12: 'e_antineutrino', 14: 'mu_neutrino', -14: 'mu_antineutrino', -2212: 'antiproton',
             -2112: 'antineutron', -311: 'antikaon0', 221: 'eta', -221: 'antieta', 331: "eta'", 3122: 'lambda0',
             -3122: 'antilambda0', 3222: 'sigma+', -3222: 'antisigma+', 3212: 'sigma0', -3212: 'antisigma0', 3112: 'sigma-',
             -3112: 'antisigma-', 3322: 'xi0', -3322: 'antixi0', 3312: 'xi-', -3312: 'antixi-', 3334: 'omega-',
             -3334: 'antiomega-'}

part_rev = {v: k for k, v in particles.items()}
elsymbol_rev = {v[0]: k for k, v in elements.items()}
elname_rev = {v[1]: k for k, v in elements.items()}

# Read in an ASCII dump file produced by a PHITS tally
def kf_decode(n: int) -> str:
    """Given a kf-code of a particle, return a human-readable string description."""

    if n in particles:
        return particles[n]
    elif n > 1000000:
        a = int(str(n)[-5:])
        z = (n - a) / 1000000
        return f"{a}{elements[z]}"
    else:
        raise ValueError(f"Invalid kf-code {n}.")

def kf_encode(part: str) -> int:
    """Given a particle name, return the kf-code."""
    assert isinstance(part, str), f"Invalid particle {part}. 1"
    if part in part_rev:
        return part_rev[part]
    elif len(re.split("-", part)) == 2: # element-weight
        pts = re.split("-", part)
        assert len(pts) == 2, f"Invalid particle {part}. 2"
        if pts[0] in elsymbol_rev:
            elt = elsymbol_rev[pts[0]]
            weight = int(pts[1])
            return elt * 1_000_000 + weight
        elif pts[0] in elname_rev:
            elt = elname_rev[pts[0]]
            weight = int(pts[1])
            return elt * 1_000_000 + weight
        else:
            raise ValueError(f"Invalid particle {part}. 3")
    elif m := re.match(r"([1-9][0-9]{,2})([a-zA-Z]+)", part):
        if m[2] in elsymbol_rev:
            elt = elsymbol_rev[m[2]]
            weight = int(m[1])
            return elt * 1_000_000 + weight
        elif m[2] in elname_rev:
            elt = elname_rev[m[2]]
            weight = int(m[1])
        else:
            raise ValueError(f"Invalid particle {part}. 4")
    else:
        raise ValueError(f"Invalid particle {part}. 5")




class Particle(ValSpec):
    def __init__(self):
        super().__init__(sampled_from(list(particles.values())))

    def phits(self, val):
        try:
            assert val in particles.values(), f"Invalid particle {val}."
            return kf_encode(val)
        except (AssertionError, ValueError):
            return partial(lambda va, var: ValueError(f"`{var}` must be a valid particle; got {va}"), val)

    def python(self, val):
        return val

    def description(self):
        return "a particle name"



# TODO: lack of  validity checking here-down
# class Element(ValSpec):
#     def __init__(self):
#         super().__init__(sampled_from(list(map(lambda x: x[1][0], elements.items()))))

#     def phits(self, val):
#         return val

#     def python(self, val):
#         return val

#     def description(self):
#         return "any valid element symbol/name"


class Nuclide(ValSpec):
    def __init__(self):
        @composite
        def symbol_hyphen_weight(draw):
            sym = draw(sampled_from(list(map(lambda x: x[1][0], elements.items()))))
            weight = draw(integers(min_value=1, max_value=294))
            return f"{sym}-{weight}"

        @composite
        def weight_then_symbol(draw):
            sym = draw(sampled_from(list(map(lambda x: x[1][0], elements.items()))))
            weight = draw(integers(min_value=1, max_value=294))
            return f"{weight}{sym}"

        super().__init__(one_of(symbol_hyphen_weight(), weight_then_symbol())) # in principle, this'd be cleaner were it a class attribute.

    def phits(self, val):
        try:
            assert val not in part_rev, f"Invalid nuclide {val}."
            return kf_encode(val)
        except (AssertionError, ValueError):
            return partial(lambda va, var: ValueError(f"`{var}` must be a valid nuclide; got {va}"), val)


    def python(self, val):
        return val

    def description(self):
        return "a nucleide in the form 208Pb, 208Lead, Pb-208, or Lead-208"

chemicals = ["H20", "CO2", "NH2", "NH3", "SF6", "TeF6", "CH4", "CH3", "C2H2", "C2H4", "C2H6", "C6H6", "CH32N3"]
class Chemical(ValSpec):
    def __init__(self):
        super().__init__(sampled_from(chemicals))

    def phits(self, val):
        if val in chemicals:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a valid chemical; got {va}."), val)

    def python(self, val):
        return val

    def description(self):
        return f"one of the chemical symbols {chemicals}"

class Orientation(ValSpec):
    def __init__(self):
        super().__init__(one_of(just("<"), just(">")))

    def phits(self, val):
        if val in ["<", ">"]:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be either `<` or `>`; got {va}."), val)

class Path(ValSpec):
    def __init__(self):
        super().__init__(text(alphabet=characters(whitelist_categories=("Nd", "L"))))

    def phits(self, val):
        if isinstance(val, str):
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a filename string; got {val}."), val)

    def description(self):
        return "a valid file name"

class OrthogonalMatrix(ValSpec):
    def __init__(self):
        listify = lambda x: list(map(tuple, ortho_group.rvs(x)))
        super().__init__(builds(listify, just(3)))

    def phits(self, val):
        if isinstance(val, np.ndarray):
            val = list(map(tuple, val))

        rounding = np.finfo(np.float32).eps
        if (np.array(val) @ np.transpose(val) - np.identity(3) < rounding).all():
            return val
        else:
            return partial(lambda va, var: ValueError(f"{var} must be a valid rotation matrix; got {va}."), val)

    def description(self):
        return "an orthogonal matrix (AA^T = I) representing a rotation"


class Text(ValSpec):
    def __init__(self):
        super().__init__(text())

    def phits(self, val):
        if isinstance(val, str):
            return "{" + val + "}" # TODO: un-specialize if unnecessary
        else:
            return partial(lambda va, var: ValueError(f"`{va}` must be a string; got {val}."), val)

    def description(self):
        return "a string"
