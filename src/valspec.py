from functools import partial
from hypothesis.strategies import *
from hypothesis import assume
import re
from scipy.stats import ortho_group
import numpy as np
from numpy.linalg import det

# TODO: think about if the python() methods are necessary



class ValSpec():
    # in principle, this'd be cleaner were it a class attribute.
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
        super().__init__(floats(allow_nan=False, allow_infinity=False, allow_subnormal=False, width=16))
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
        super().__init__(floats(min_value=0, exclude_min=True, allow_nan=False, allow_infinity=False, allow_subnormal=False,
                                width=16))
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
                                                                      allow_subnormal=False, width=16)))
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
        super().__init__(floats(min_value=start, max_value=stop, allow_nan=False, allow_infinity=False, allow_subnormal=False,
                                width=16))
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


@composite
def builds_right(draw, cl, re, op):
    try:
        ob = draw(builds(cl, *re, **op))
    except ValueError as e:
        ob = None

    assume(ob)
    return ob

class IsA(ValSpec):
    def __init__(self, cls, index=False):
        req = []
        for phits_iden, valspec, idx, *s in sorted((v for v in cls.syntax.values() if v[2] is not None), key=lambda t: t[2]):
            if isinstance(valspec, tuple):
                req.append(tuples(*[i.strat for i in valspec]))
            else:
                req.append(valspec.strat)

        opt = dict()
        for py_iden, (phits_iden, valspec, idx, *s) in filter(lambda t: t[1][2] is None, cls.syntax.items()):
            if isinstance(valspec, tuple):
                opt[py_iden] = one_of(none(), tuples(*[i.strat for i in valspec]))
            else:
                opt[py_iden] = one_of(none(), valspec.strat)


        super().__init__(builds_right(cls, req, opt)) # TODO: sus
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



# Generated from the file $PHITSHOME/data/xsdir.jnd; <atomic number>: ('symbol', 'name', ((<atomic_weight>, 'library_code'))).
# Atomic weight of zero I presume corresponds to a measurement of an average over atomic weights.
elements = {1: ('H', 'Hydrogen', ((1, '50c'), (2, '50c'), (0, '50p'), (0, '50e'), (1, '51c'), (2, '51c'), (1, '51h'), (2, '51h'))),
            2: ('He', 'Helium', ((3, '50c'), (4, '50c'), (0, '50p'), (0, '50e'))),
            3: ('Li', 'Lithium', ((6, '50c'), (7, '50c'), (0, '50p'), (0, '50e'), (6, '51h'), (7, '51h'))),
            4: ('Be', 'Beryllium', ((9, '50c'), (0, '50p'), (0, '50e'), (9, '51h'))),
            5: ('B', 'Boron', ((10, '50c'), (11, '50c'), (0, '50p'), (0, '50e'))),
            6: ('C', 'Carbon', ((0, '50c'), (12, '50c'), (0, '50p'), (0, '50e'), (0, '51c'), (12, '51c'), (13, '51c'), (0, '51h'),
                                (12, '51h'), (13, '51h'))),
            7: ('N', 'Nitrogen', ((14, '50c'), (15, '50c'), (0, '50p'), (0, '50e'), (14, '51c'), (14, '51h'))),
            8: ('O', 'Oxygen', ((16, '50c'), (0, '50p'), (0, '50e'), (16, '51c'), (16, '51h'))),
            9: ('F', 'Fluorine', ((19, '50c'), (0, '50p'), (0, '50e'))),
            10: ('Ne', 'Neon', ((0, '50p'), (0, '50e'))),
            11: ('Na', 'Sodium', ((23, '50c'), (0, '50p'), (0, '50e'))),
            12: ('Mg', 'Magnesium', ((24, '50c'), (25, '50c'), (26, '50c'), (0, '50p'), (0, '50e'))),
            13: ('Al', 'Aluminum', ((27, '50c'), (0, '50p'), (0, '50e'), (27, '51c'), (27, '51h'))),
            14: ('Si', 'Silicon', ((28, '50c'), (29, '50c'), (30, '50c'), (0, '50p'), (0, '50e'), (28, '51c'), (29, '51c'), (30, '51c'),
                                   (28, '51h'), (29, '51h'), (30, '51h'))),
            15: ('P', 'Phosphorus', ((31, '50c'), (0, '50p'), (0, '50e'))),
            16: ('S', 'Sulfur', ((32, '50c'), (33, '50c'), (34, '50c'), (36, '50c'), (0, '50p'), (0, '50e'))),
            17: ('Cl', 'Chlorine', ((35, '50c'), (37, '50c'), (0, '50p'), (0, '50e'))),
            18: ('Ar', 'Argon', ((40, '50c'), (0, '50p'), (0, '50e'))),
            19: ('K', 'Potassium', ((39, '50c'), (40, '50c'), (41, '50c'), (0, '50p'), (0, '50e'))),
            20: ('Ca', 'Calcium', ((40, '50c'), (42, '50c'), (43, '50c'), (44, '50c'), (46, '50c'), (48, '50c'), (0, '50p'), (0, '50e'))),
            21: ('Sc', 'Scandium', ((45, '50c'), (0, '50p'), (0, '50e'))),
            22: ('Ti', 'Titanium', ((46, '50c'), (47, '50c'), (48, '50c'), (49, '50c'), (50, '50c'), (0, '50p'), (0, '50e'))),
            23: ('V', 'Vanadium', ((50, '50c'), (51, '50c'), (0, '50p'), (0, '50e'))),
            24: ('Cr', 'Chromium', ((50, '50c'), (52, '50c'), (53, '50c'), (54, '50c'), (0, '50p'), (0, '50e'))),
            25: ('Mn', 'Manganese', ((55, '50c'), (0, '50p'), (0, '50e'))),
            26: ('Fe', 'Iron', ((54, '50c'), (56, '50c'), (57, '50c'), (58, '50c'), (59, '50c'), (0, '50p'), (0, '50e'), (54, '51c'),
                                (56, '51c'), (57, '51c'), (58, '51c'), (54, '51h'), (56, '51h'), (57, '51h'), (58, '51h'))),
            27: ('Co', 'Cobalt', ((59, '50c'), (0, '50p'), (0, '50e'))),
            28: ('Ni', 'Nickel', ((58, '50c'), (59, '50c'), (60, '50c'), (61, '50c'), (62, '50c'), (64, '50c'), (0, '50p'), (0, '50e'))),
            29: ('Cu', 'Copper', ((63, '50c'), (65, '50c'), (0, '50p'), (0, '50e'), (63, '51c'), (65, '51c'), (63, '51h'), (65, '51h'))),
            30: ('Zn', 'Zinc', ((64, '50c'), (65, '50c'), (66, '50c'), (67, '50c'), (68, '50c'), (70, '50c'), (0, '50p'), (0, '50e'))),
            31: ('Ga', 'Gallium', ((69, '50c'), (71, '50c'), (0, '50p'), (0, '50e'))),
            32: ('Ge', 'Germanium', ((70, '50c'), (72, '50c'), (73, '50c'), (74, '50c'), (76, '50c'), (0, '50p'), (0, '50e'))),
            33: ('As', 'Arsenic', ((75, '50c'), (0, '50p'), (0, '50e'))),
            34: ('Se', 'Selenium', ((74, '50c'), (76, '50c'), (77, '50c'), (78, '50c'), (79, '50c'), (80, '50c'), (82, '50c'), (0, '50p'),
                                    (0, '50e'))),
            35: ('Br', 'Bromine', ((79, '50c'), (81, '50c'), (0, '50p'), (0, '50e'))),
            36: ('Kr', 'Krypton', ((78, '50c'), (80, '50c'), (82, '50c'), (83, '50c'), (84, '50c'), (85, '50c'), (86, '50c'), (0, '50p'),
                                   (0, '50e'))),
            37: ('Rb', 'Rubidium', ((85, '50c'), (86, '50c'), (87, '50c'), (0, '50p'), (0, '50e'))),
            38: ('Sr', 'Strontium', ((84, '50c'), (86, '50c'), (87, '50c'), (88, '50c'), (89, '50c'), (90, '50c'), (0, '50p'),
                                     (0, '50e'))),
            39: ('Y', 'Yttrium', ((89, '50c'), (90, '50c'), (91, '50c'), (0, '50p'), (0, '50e'))),
            40: ('Zr', 'Zirconium', ((90, '50c'), (91, '50c'), (92, '50c'), (93, '50c'), (94, '50c'), (95, '50c'), (96, '50c'), (0, '50p'),
                                     (0, '50e'))),
            41: ('Nb', 'Niobium', ((93, '50c'), (94, '50c'), (95, '50c'), (0, '50p'), (0, '50e'))),
            42: ('Mo', 'Molybdenum', ((92, '50c'), (94, '50c'), (95, '50c'), (96, '50c'), (97, '50c'), (98, '50c'), (99, '50c'),
                                      (100, '50c'), (0, '50p'), (0, '50e'))),
            43: ('Tc', 'Technetium', ((99, '50c'), (0, '50p'), (0, '50e'))),
            44: ('Ru', 'Ruthenium', ((96, '50c'), (98, '50c'), (99, '50c'), (100, '50c'), (101, '50c'), (102, '50c'), (103, '50c'),
                                     (104, '50c'), (105, '50c'), (106, '50c'), (0, '50p'), (0, '50e'))),
            45: ('Rh', 'Rhodium', ((103, '50c'), (105, '50c'), (0, '50p'), (0, '50e'))),
            46: ('Pd', 'Palladium', ((102, '50c'), (104, '50c'), (105, '50c'), (106, '50c'), (107, '50c'), (108, '50c'), (110, '50c'),
                                     (0, '50p'), (0, '50e'))),
            47: ('Ag', 'Silver', ((107, '50c'), (109, '50c'), (190, '50c'), (111, '50c'), (0, '50p'), (0, '50e'))),
            48: ('Cd', 'Cadmium', ((106, '50c'), (108, '50c'), (110, '50c'), (111, '50c'), (112, '50c'), (113, '50c'), (114, '50c'),
                                   (116, '50c'), (0, '50p'), (0, '50e'))),
            49: ('In', 'Indium', ((113, '50c'), (115, '50c'), (0, '50p'), (0, '50e'))),
            50: ('Sn', 'Tin', ((112, '50c'), (114, '50c'), (115, '50c'), (116, '50c'), (117, '50c'), (118, '50c'), (119, '50c'),
                               (120, '50c'), (122, '50c'), (123, '50c'), (124, '50c'), (126, '50c'), (0, '50p'), (0, '50e'))),
            51: ('Sb', 'Antimony', ((121, '50c'), (123, '50c'), (124, '50c'), (125, '50c'), (126, '50c'), (0, '50p'), (0, '50e'))),
            52: ('Te', 'Tellurium', ((120, '50c'), (122, '50c'), (123, '50c'), (124, '50c'), (125, '50c'), (126, '50c'), (197, '50c'),
                                     (128, '50c'), (199, '50c'), (130, '50c'), (132, '50c'), (0, '50p'), (0, '50e'))),
            53: ('I', 'Iodine', ((127, '50c'), (129, '50c'), (130, '50c'), (131, '50c'), (135, '50c'), (0, '50p'), (0, '50e'))),
            54: ('Xe', 'Xenon', ((124, '50c'), (126, '50c'), (128, '50c'), (129, '50c'), (130, '50c'), (131, '50c'), (132, '50c'),
                                 (133, '50c'), (134, '50c'), (135, '50c'), (136, '50c'), (0, '50p'), (0, '50e'))),
            55: ('Cs', 'Cesium', ((133, '50c'), (134, '50c'), (135, '50c'), (136, '50c'), (137, '50c'), (0, '50p'), (0, '50e'))),
            56: ('Ba', 'Barium', ((130, '50c'), (132, '50c'), (134, '50c'), (135, '50c'), (136, '50c'), (137, '50c'), (138, '50c'),
                                  (140, '50c'), (0, '50p'), (0, '50e'))),
            57: ('La', 'Lanthanum', ((138, '50c'), (139, '50c'), (140, '50c'), (0, '50p'), (0, '50e'))),
            58: ('Ce', 'Cerium', ((140, '50c'), (141, '50c'), (142, '50c'), (143, '50c'), (144, '50c'), (0, '50p'), (0, '50e'))),
            59: ('Pr', 'Praseodymium', ((141, '50c'), (143, '50c'), (0, '50p'), (0, '50e'))),
            60: ('Nd', 'Neodymium', ((142, '50c'), (143, '50c'), (144, '50c'), (145, '50c'), (146, '50c'), (147, '50c'), (148, '50c'),
                                     (150, '50c'), (0, '50p'), (0, '50e'))),
            61: ('Pm', 'Promethium', ((147, '50c'), (148, '50c'), (198, '50c'), (149, '50c'), (151, '50c'), (0, '50p'), (0, '50e'))),
            62: ('Sm', 'Samarium', ((144, '50c'), (147, '50c'), (148, '50c'), (149, '50c'), (150, '50c'), (151, '50c'), (152, '50c'),
                                    (153, '50c'), (154, '50c'), (0, '50p'), (0, '50e'))),
            63: ('Eu', 'Europium', ((151, '50c'), (152, '50c'), (153, '50c'), (154, '50c'), (155, '50c'), (156, '50c'), (157, '50c'),
                                    (0, '50p'), (0, '50e'))),
            64: ('Gd', 'Gadolinium', ((152, '50c'), (153, '50c'), (154, '50c'), (155, '50c'), (156, '50c'), (157, '50c'), (158, '50c'),
                                      (160, '50c'), (0, '50p'), (0, '50e'))),
            65: ('Tb', 'Terbium', ((159, '50c'), (160, '50c'), (0, '50p'), (0, '50e'))),
            66: ('Dy', 'Dysprosium', ((154, '50c'), (156, '50c'), (158, '50c'), (159, '50c'), (160, '50c'), (161, '50c'), (162, '50c'),
                                      (163, '50c'), (164, '50c'), (0, '50p'), (0, '50e'))),
            67: ('Ho', 'Holmium', ((0, '50p'), (0, '50e'))),
            68: ('Er', 'Erbium', ((162, '50c'), (164, '50c'), (166, '50c'), (167, '50c'), (168, '50c'), (170, '50c'), (0, '50p'),
                                  (0, '50e'))),
            69: ('Tm', 'Thulium', ((169, '50c'), (0, '50p'), (0, '50e'))),
            70: ('Yb', 'Ytterbium', ((168, '50c'), (170, '50c'), (171, '50c'), (172, '50c'), (173, '50c'), (174, '50c'), (176, '50c'),
                                     (0, '50p'), (0, '50e'))),
            71: ('Lu', 'Lutetium', ((0, '50p'), (0, '50e'))),
            72: ('Hf', 'Hafnium', ((174, '50c'), (176, '50c'), (177, '50c'), (178, '50c'), (179, '50c'), (180, '50c'), (181, '50c'),
                                   (182, '50c'), (0, '50p'), (0, '50e'))),
            73: ('Ta', 'Tantalum', ((181, '50c'), (0, '50p'), (0, '50e'))),
            74: ('W', 'Tungsten', ((180, '50c'), (182, '50c'), (183, '50c'), (184, '50c'), (186, '50c'), (0, '50p'), (0, '50e'))),
            75: ('Re', 'Rhenium', ((0, '50p'), (0, '50e'))),
            76: ('Os', 'Osmium', ((184, '50c'), (186, '50c'), (187, '50c'), (188, '50c'), (189, '50c'), (190, '50c'), (192, '50c'),
                                  (0, '50p'), (0, '50e'))),
            77: ('Ir', 'Iridium', ((0, '50p'), (0, '50e'))),
            78: ('Pt', 'Platinum', ((0, '50p'), (0, '50e'))),
            79: ('Au', 'Gold', ((197, '50c'), (0, '50p'), (0, '50e'))),
            80: ('Hg', 'Mercury', ((196, '50c'), (198, '50c'), (199, '50c'), (200, '50c'), (201, '50c'), (202, '50c'), (204, '50c'),
                                   (0, '50p'), (0, '50e'))),
            81: ('Tl', 'Thallium', ((0, '50p'), (0, '50e'))),
            82: ('Pb', 'Lead', ((204, '50c'), (206, '50c'), (207, '50c'), (208, '50c'), (0, '50p'), (0, '50e'), (204, '51c'), (206, '51c'),
                                (207, '51c'), (208, '51c'), (204, '51h'), (206, '51h'), (207, '51h'), (208, '51h'))),
            83: ('Bi', 'Bismuth', ((209, '50c'), (0, '50p'), (0, '50e'), (209, '51c'), (209, '51h'))),
            84: ('Po', 'Polonium', ((0, '50p'), (0, '50e'))),
            85: ('At', 'Astatine', ((0, '50p'), (0, '50e'))),
            86: ('Rn', 'Radon', ((0, '50p'), (0, '50e'))),
            87: ('Fr', 'Francium', ((0, '50p'), (0, '50e'))),
            88: ('Ra', 'Radium', ((223, '50c'), (224, '50c'), (225, '50c'), (226, '50c'), (0, '50p'), (0, '50e'))),
            89: ('Ac', 'Actinium', ((225, '50c'), (226, '50c'), (227, '50c'), (0, '50p'), (0, '50e'))),
            90: ('Th', 'Thorium', ((227, '50c'), (228, '50c'), (229, '50c'), (230, '50c'), (231, '50c'), (232, '50c'), (233, '50c'),
                                   (234, '50c'), (0, '50p'), (0, '50e'))),
            91: ('Pa', 'Protactinium', ((229, '50c'), (230, '50c'), (231, '50c'), (232, '50c'), (233, '50c'), (0, '50p'), (0, '50e'))),
            92: ('U', 'Uranium', ((230, '50c'), (231, '50c'), (232, '50c'), (233, '50c'), (234, '50c'), (235, '50c'), (236, '50c'),
                                  (237, '50c'), (238, '50c'), (0, '50p'), (0, '50e'))),
            93: ('Np', 'Neptunium', ((234, '50c'), (235, '50c'), (236, '50c'), (237, '50c'), (238, '50c'), (239, '50c'), (0, '50p'),
                                     (0, '50e'))),
            94: ('Pu', 'Plutonium', ((236, '50c'), (237, '50c'), (238, '50c'), (239, '50c'), (240, '50c'), (241, '50c'), (242, '50c'),
                                     (244, '50c'), (246, '50c'), (0, '50p'), (0, '50e'))),
            95: ('Am', 'Americium', ((240, '50c'), (241, '50c'), (242, '50c'), (292, '50c'), (243, '50c'), (244, '50c'), (294, '50c'),
                                     (0, '50p'), (0, '50e'))),
            96: ('Cm', 'Curium', ((240, '50c'), (241, '50c'), (242, '50c'), (243, '50c'), (244, '50c'), (245, '50c'), (246, '50c'),
                                  (247, '50c'), (248, '50c'), (249, '50c'), (250, '50c'), (0, '50p'), (0, '50e'))),
            97: ('Bk', 'Berkelium', ((245, '50c'), (246, '50c'), (247, '50c'), (248, '50c'), (249, '50c'), (250, '50c'), (0, '50p'),
                                     (0, '50e'))),
            98: ('Cf', 'Californium', ((246, '50c'), (248, '50c'), (249, '50c'), (250, '50c'), (251, '50c'), (252, '50c'), (253, '50c'),
                                       (254, '50c'), (0, '50p'), (0, '50e'))),
            99: ('Es', 'Einsteinium', ((251, '50c'), (252, '50c'), (253, '50c'), (254, '50c'), (294, '50c'), (255, '50c'), (0, '50p'),
                                       (0, '50e'))),
            100: ('Fm', 'Fermium', ((255, '50c'), (0, '50p'), (0, '50e')))}

particles = {2212: "proton", 2112: "neutron", 211: "pion+", 111: "pion0", -211: "pion-", -13: "muon+", 13: "muon-",
             321: "kaon+", 311: "kaon0", -321: "kaon-", 11: "electron", -11: "positron", 22: "photon",
             12: 'other', -12: 'other', 14: 'other', -14: 'other', -2212: 'other',
             -2112: 'other', -311: 'other', 221: 'other', -221: 'other', 331: 'other', 3122: 'other',
             -3122: 'other', 3222: 'other', -3222: 'other', 3212: 'other', -3212: 'other', 3112: 'other',
             -3112: 'other', 3322: 'other', -3322: 'other', 3312: 'other', -3312: 'other', 3334: 'other',
             -3334: 'other'}

part_rev = {v: k for k, v in particles.items()}
elsymbol_rev = {v[0]: k for k, v in elements.items()}
elname_rev = {v[1]: k for k, v in elements.items()}
el_weights = {k: set(map(lambda x: x[0], v[2])) for k, v in elements.items()}
j4_el_weights = {k: set(map(lambda x: x[0], v[2])) for k, v in \
                 filter(lambda x: set(("50c", "51h")).issubset(set(map(lambda x: x[1], x[1][2]))),
                        elements.items())}

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
    assert isinstance(part, str), f"Invalid particle {part}; must be a string."
    if part in part_rev:
        return part_rev[part]

    elif part in elsymbol_rev:
        return int(elsymbol_rev[part]) * 1_000

    elif part in elname_rev:
        return int(elname_rev[part]) * 1_000

    elif len(re.split("-", part)) == 2: # element-weight
        pts = re.split("-", part)
        assert len(pts) == 2, f"Invalid particle {part}; must be of the form 208Pb or Pb-208 (more than two components)."
        capped = pts[0].title()
        if capped in elsymbol_rev:
            elt = elsymbol_rev[capped]
            weight = int(pts[1])
            assert weight in el_weights[elt], f"Unsupported isotope {part}; acceptable values are {el_weights[elt]}."
            return f"{elt}{weight:03d}"
        elif capped in elname_rev:
            elt = elname_rev[capped]
            weight = int(pts[1])
            assert weight in el_weights[elt], f"Unsupported isotope {part}; acceptable atomic weights are {el_weights[elt]}."
            return f"{elt}{weight:03d}"
        else:
            raise ValueError(f"Invalid particle {part}; must be of the form 208Pb or Pb-208 (symbol or name not found).")

    elif m := re.match(r"([1-9][0-9]{,2})([a-zA-Z]+)", part):
        capped = m[2].title()
        if capped in elsymbol_rev:
            elt = elsymbol_rev[capped]
            weight = int(m[1])
            return f"{elt}{weight:03d}"
        elif capped in elname_rev:
            elt = elname_rev[capped]
            weight = int(m[1])
            return f"{elt}{weight:03d}"
        else:
            raise ValueError(f"Invalid particle {part}; must be of the form 208Pb or Pb-208 (symbol or name not found).")

    else:
        raise ValueError(f"Invalid particle {part}.")




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
            wgts = list(filter(lambda x: x != 0,
                               map(lambda x: x[0], elements[elsymbol_rev[sym]][2])))
            assume(len(wgts) > 0)
            weight = draw(sampled_from(list(wgts)))
            return f"{sym}-{weight}"

        @composite
        def weight_then_symbol(draw):
            sym = draw(sampled_from(list(map(lambda x: x[1][0], elements.items()))))
            wgts = list(filter(lambda x: x != 0,
                               map(lambda x: x[0], elements[elsymbol_rev[sym]][2])))
            assume(len(wgts) > 0)
            weight = draw(sampled_from(list(wgts)))
            return f"{weight}{sym}"

        @composite
        def symbol_alone(draw):
            sym = draw(sampled_from(list(map(lambda x: x[1][0], elements.items()))))
            assume(0 in map(lambda x: x[0], elements[elsymbol_rev[sym]][2]))
            return f"{sym}"

        @composite
        def name_hyphen_weight(draw):
            name = draw(sampled_from(list(map(lambda x: x[1][1], elements.items()))))
            wgts = list(filter(lambda x: x != 0,
                               map(lambda x: x[0], elements[elname_rev[name]][2])))
            assume(len(wgts) > 0)
            weight = draw(sampled_from(list(wgts)))
            return f"{name}-{weight}"

        @composite
        def weight_then_name(draw):
            name = draw(sampled_from(list(map(lambda x: x[1][1], elements.items()))))

            wgts = list(filter(lambda x: x != 0,
                               map(lambda x: x[0], elements[elname_rev[name]][2])))
            assume(len(wgts) > 0)
            weight = draw(sampled_from(list(wgts)))
            return f"{weight}{name}"

        @composite
        def name_alone(draw):
            name = draw(sampled_from(list(map(lambda x: x[1][1], elements.items()))))
            assume(0 in map(lambda x: x[0], elements[elname_rev[name]][2]))
            return f"{name}"




        super().__init__(one_of(symbol_hyphen_weight(), weight_then_symbol(), symbol_alone(),

                                name_hyphen_weight(), weight_then_name(), name_alone()))

    def phits(self, val):
        try:
            assert val not in part_rev, f"Particle {val} is not a nuclide."
            return kf_encode(val)
        except (AssertionError, ValueError) as e:

            return partial(lambda va, er, var: ValueError(f"`{var}` must be a valid nuclide; {va} resulted in error:\n{er}"),
                           val, e)


    def python(self, val):
        return val

    def description(self):
        return "a nucleide in the form 208Pb, 208Lead, Pb-208, or Lead-208"


class JENDL4Nuclide(ValSpec):
    def __init__(self):
        @composite
        def symbol_hyphen_weight(draw):
            sym = draw(sampled_from(list(map(lambda x: x[1][0], elements.items()))))
            wgts = list(filter(lambda x: x != 0,
                               map(lambda x: x[0], elements[elsymbol_rev[sym]][2])))
            assume(len(wgts) > 0) # TODO: necessary?
            weight = draw(sampled_from(list(wgts)))
            return f"{sym}-{weight}"

        @composite
        def weight_then_symbol(draw):
            sym = draw(sampled_from(list(map(lambda x: x[1][0], elements.items()))))
            wgts = list(filter(lambda x: x != 0,
                               map(lambda x: x[0], elements[elsymbol_rev[sym]][2])))
            assume(len(wgts) > 0)
            weight = draw(sampled_from(list(wgts)))
            return f"{weight}{sym}"

        @composite
        def symbol_alone(draw):
            sym = draw(sampled_from(list(map(lambda x: x[1][0], elements.items()))))
            assume(elsymbol_rev[sym] in j4_el_weights)
            assume(0 in map(lambda x: x[0], elements[elsymbol_rev[sym]][2])) # TODO: necessary?
            return f"{sym}"

        @composite
        def name_hyphen_weight(draw):
            name = draw(sampled_from(list(map(lambda x: x[1][1], elements.items()))))
            wgts = list(filter(lambda x: x != 0,
                               map(lambda x: x[0], elements[elname_rev[name]][2])))
            assume(len(wgts) > 0)
            weight = draw(sampled_from(list(wgts)))
            return f"{name}-{weight}"

        @composite
        def weight_then_name(draw):
            name = draw(sampled_from(list(map(lambda x: x[1][1], elements.items()))))

            wgts = list(filter(lambda x: x != 0,
                               map(lambda x: x[0], elements[elname_rev[name]][2])))
            assume(len(wgts) > 0)
            weight = draw(sampled_from(list(wgts)))
            return f"{weight}{name}"

        @composite
        def name_alone(draw):
            name = draw(sampled_from(list(map(lambda x: x[1][1], elements.items()))))
            assume(elname_rev[name] in j4_el_weights)
            assume(0 in map(lambda x: x[0], elements[elname_rev[name]][2])) # TODO: necessary?
            return f"{name}"




        super().__init__(one_of(symbol_hyphen_weight(), weight_then_symbol(), symbol_alone(),

                                name_hyphen_weight(), weight_then_name(), name_alone()))

    def phits(self, val):
        try:
            assert val not in part_rev, f"Particle {val} is not a nuclide."
            return kf_encode(val)
        except (AssertionError, ValueError) as e:

            return partial(lambda va, er, var: ValueError(f"`{var}` must be a valid nuclide; {va} resulted in error:\n{er}"),
                           val, e)


    def python(self, val):
        return val

    def description(self):
        return "a nucleide in the form 208Pb, 208Lead, Pb-208, or Lead-208"


material_libs = ['lmeth.20t', 'lwtr.20t', 'be.24t', 'smeth.20t', 'grph.22t', 'grph.26t', 'beo.25t', 'benz.25t', 'beo.21t', 'hwtr.25t',
                 'grph.29t', 'beo.24t', 'be.20t', 'zr_h.25t', 'hwtr.21t', 'dpara.20t', 'grph.28t', 'zr_h.27t', 'be.22t', 'beo.23t',
                 'benz.27t', 'grph.25t', 'lwtr.22t', 'grph.27t', 'hwtr.24t', 'h_zr.26t', 'h_zr.23t', 'h_zr.20t', 'grph.23t', 'be.27t',
                 'poly.21t', 'benz.26t', 'hwtr.23t', 'grph.20t', 'lwtr.21t', 'grph.21t', 'zr_h.22t', 'h_zr.24t', 'benz.24t', 'hwtr.27t',
                 'beo.22t', 'hwtr.22t', 'zr_h.24t', 'be.26t', 'lwtr.25t', 'be.25t', 'benz.20t', 'dortho.20t', 'hortho.20t', 'benz.21t',
                 'lwtr.23t', 'beo.26t', 'lwtr.27t', 'hwtr.20t', 'grph.24t', 'beo.20t', 'zr_h.20t', 'be.21t', 'beo.27t', 'h_zr.22t',
                 'zr_h.21t', 'lwtr.26t', 'benz.22t', 'h_zr.25t', 'benz.23t', 'zr_h.26t', 'hpara.20t', 'be.23t', 'hwtr.26t', 'lwtr.24t',
                 'poly.20t', 'zr_h.23t', 'h_zr.27t', 'h_zr.21t']

class ThermalLib(ValSpec):
    def __init__(self):
        super().__init__(sampled_from(material_libs))

    def phits(self, val):
        if val in material_libs:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a valid material library; got {va}."), val)

    def description(self):
        return "a material thermal neutron scattering law library identifier"

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
        super().__init__(text(min_size=1, alphabet=characters(min_codepoint=0x0041, max_codepoint=0x007a)))

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

        rounding = 2.0e-6
        if (np.array(val) @ np.transpose(val) - np.identity(3) < rounding).all():
            return val
        else:
            return partial(lambda va, var: ValueError(f"{var} must be a valid rotation matrix; got {va}."), val)

    def description(self):
        return "an orthogonal matrix (AA^T = I) representing a rotation"


class LibraryID(ValSpec):
    def __init__(self):
        # TODO: make sure I understand the libraries right
        self.ids = ["50c", "20t", "21t", "22t", "23t", "24t", "25t", "26t", "27t", "28t", "29t", "50t", "50e", "51c", "51h", "90p"]
        super().__init__(sampled_from(self.ids))

    def phits(self, val):
        if val in self.ids:
            return val
        else:
            return partial(lambda va, var: ValueError(f"{var} must be a valid library ID (one of {self.ids}); got {va}."), val)

    def description(self):
        return f"a library ID (one of {self.ids})"

class Text(ValSpec):
    def __init__(self):

        super().__init__(text(min_size=1, alphabet=characters(min_codepoint=0x0041, max_codepoint=0x007a))) # TODO: make more general

    def phits(self, val):
        if isinstance(val, str) and val != "":
            return "{ " + val + " }" # TODO: un-specialize if unnecessary
        else:
            return partial(lambda va, var: ValueError(f"`{va}` must be a string; got {val}."), val)

    def description(self):
        return "a string"


named_colors = ["white", "lightgray", "gray", "darkgray", "matblack", "black", "red", "orange", "yellow", "green", "cyan", "blue",
                "violet", "magenta", "darkred", "pink", "pastelpink", "orange", "brown", "darkbrown", "pastelbrown", "orangeyellow",
                "camel", "pastelyellow", "yellow", "pastelgreen", "yellowgreen", "green", "darkgreen", "mossgreen", "bluegreen",
                "pastelcyan", "pastelblue", "cyan", "cyanblue", "violet", "purple", "magenta", "winered", "pastelmagenta",
                "pastelpurple", "pastelviolet"]
class Color(ValSpec):
    def __init__(self):
        super().__init__(one_of(sampled_from(named_colors),
                                tuples(floats(min_value=0, max_value=1, allow_nan=False, allow_infinity=False, allow_subnormal=False),
                                       floats(min_value=0, max_value=1, allow_nan=False, allow_infinity=False, allow_subnormal=False),
                                       floats(min_value=0, max_value=1, allow_nan=False, allow_infinity=False, allow_subnormal=False))))

    def phits(self, val):
        if val in named_colors :
            return val
        elif isinstance(val, tuple) and len(val) == 3 and all(map(lambda x: 0 <= x <= 1, val)):
            return f"{{ {val[0]} {val[1]} {val[2]} }}"
        else:
            return partial(lambda va, var: ValueError(f"`{va}` must be a valid color; got {val}."), val)

    def description(self):
        return "a color"
