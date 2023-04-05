from functools import partial
from hypothesis.strategies import *
import re
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
        return "int literal"

class Real(ValSpec):
    def __init__(self):
        super().__init__(floats())
    def phits(self, val):
        if isinstance(val, float):
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a float; got {va}."), val)

    def python(self, val):
        return val

    def description(self):
        return "float literal"

class PosInt(ValSpec):
    def __init__(self):
        super().__init__(integers(min_value=0))

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
        return "int literal strictly greater than zero"


class PosReal(ValSpec):
    def __init__(self):
        super().__init__(floats(min_value=0))
    def phits(self, val):
        if isinstance(val, float) and val > 0:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a positive float; got {va}."), val)

    def python(self, val):
        if val > 0:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be positive; got {va}."), val)

    def description(self):
        return "float literal strictly greater than zero"




class NegDisable(ValSpec):
    def __init__(self):
        super().__init__(one_of(none(), integers(), floats()))
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
        return "either None or a number literal"

class Between(ValSpec):
    def __init__(self, start, stop):
        super().__init__(integers(min_vale=start, max_value=stop))
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
        return f"integer literal between {self.start} and {self.stop}"



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
        return f"one of {list(self.dic.keys())}"


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
            return val.index
        else:
            return val.definition()

    def python(self, val):
        return val

    def description(self):
        if self.index:
            return f"an instance of {self.cls}"
        else:
            return f"an instance of {self.cls}"



class List(ValSpec):
    def __init__(self, entr):
        super().__init__(lists(entr.strat))
        self.entr = entr

    def phits(self, val):
        return list(map(lambda v: self.entr.phits(v), val))

    def python(self, val):
        return list(map(lambda v: self.entr.python(v), val))

    def description(self):
        return f"a list of values, each of which is {entr.description()}"

class Tuple(ValSpec):
    def __init__(self, *entr):
        super().__init__(tuples(*[i.strat for i in entr]))
        self.entr = entr

    def phits(self, val):
        return list(map(lambda v: self.entr.phits(v), val))

    def python(self, val):
        return list(map(lambda v: self.entr.python(v), val))

    def description(self):
        return f"a list of values, each of which is {entr.description()}"


class OneOf(ValSpec):
    def __init__(self, *args):
        assert all(map(lambda x: isinstance(x, ValSpec), args)), "All arguments to OneOf must be value specifications."
        super().__init__(one_of(*[i.strat for i in args]))
        self.choices = args

    def phits(self, val):
        return self.that_which_applies(val).phits(val)

    def python(self, val):
        return self.that_which_applies(val).python(val)

    def that_which_applies(self, val):
        def _applies(s, val):
            try:
                s.phits(val)
                return ("phits", s)
            except Exception:
                try:
                    s.python(val)
                    return ("python", s)
                except Exception:
                    return None

        applicable = filter(lambda x: _applies(x, val) is not None, self.choices)
        if len(applicable) == 0:
            return ValueError("Empty OneOf value specification.")
        else:
            assert len(applicable) == 1, "Ambiguous OneOf value specification."
            return applicable[0]

    def description(self):
        return "either" + ", ".join(map(lambda x: x.description(), self.choices[:-1])) + ", or " + self.choices[-1].description()



class Particle(ValSpec):
    def __init__(self):
        self.re = r"(?i:proton|neutron|pion(+|0|-)|muon(+|-)|kaon(+|0|-)|electron|positron|photon|deuteron|triton|3he|alpha)"
        super().__init__(one_of(from_regex(self.re, fullmatch=True), Nucleon().strat, just("all")))

    def phits(self, val):
        if re.fullmatch(self.re, val):
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a valid particle name; got {va}."), val)

    def python(self, val):
        return val

el_symbols = {"ac", "ag", "al", "am", "ar", "as", "at", "au", "b", "ba", "be", "bh", "bi", "bk", "br", "c", "ca", "cd", "ce", "cf", "cl", "cm",
              "cn", "co", "cr", "cs", "cu", "db", "ds", "dy", "er", "es", "eu", "f", "fe", "fl", "fm", "fr", "ga", "gd", "ge", "h",  "he","hf",
              "hg", "ho", "hs", "i", "in", "ir", "k", "kr", "la", "li", "lr", "lu", "lv", "mc", "md", "mg", "mn", "mo", "mt", "n", "na", "nb",
              "nd", "ne", "nh", "ni", "no", "np", "o", "og", "os", "p", "pa", "pb", "pd", "pm", "po", "pr", "pt", "pu", "ra", "rb", "re", "rf",
              "rg", "rh", "rn", "ru", "s", "sb", "sc", "se", "sg", "si", "sm", "sn", "sr", "ta", "tb", "tc", "te", "th", "ti", "tl", "tm",
              "ts", "u", "v", "w", "xe", "y", "yb", "zn", "zr"}

mat_identifiers = {"actinium", "silver", "aluminum", "americium", "argon", "arsenic",
                   "astatine", "gold", "boron", "barium", "beryllium", "bohrium",
                   "bismuth", "berkelium", "bromine", "carbon", "calcium", "cadmium",
                   "cerium", "californium", "chlorine", "curium", "copernicium", "cobalt",
                   "chromium", "cesium", "copper", "dubnium", "darmstadtium", "dysprosium",
                   "erbium", "einsteinium", "europium", "fluorine", "iron", "flerovium",
                   "fermium", "francium", "gallium", "gadolinium", "germanium", "hydrogen",
                   "helium", "hafnium", "mercury", "holmium", "hassium", "iodine",
                   "indium", "iridium", "potassium", "krypton", "lanthanum", "lithium",
                   "lawrencium", "lutetium", "livermorium", "moscovium", "mendelevium",
                   "magnesium", "manganese", "molybdenum", "meitnerium", "nitrogen", "sodium",
                   "niobium", "neodymium", "neon", "nihonium", "nickel", "nobelium",
                   "neptunium", "oxygen", "oganesson", "osmium", "phosphorus", "protactinium",
                   "lead", "palladium", "promethium", "polonium", "praseodymium", "platinum",
                   "plutonium", "radium", "rubidium", "rhenium", "rutherfordium", "roentgenium", "rhodium", "radon", "ruthenium", "sulfur",
                   "antimony", "scandium", "selenium", "seaborgium", "silicon", "samarium", "tin", "strontium", "tantalum", "terbium",
                   "technetium", "tellurium", "thorium", "titanium", "thallium", "thulium", "tennessine", "uranium",
                   "vanadium", "tungsten", "xenon", "yttrium", "ytterbium", "zinc", "zirconium"}

# TODO: lack of  validity checking here-down
class Element(ValSpec):
    def __init__(self):
        super().__init__(sampled_from(el_symbols))

    def phits(self, val):
        return val

    def python(self, val):
        return val


class Nucleide(ValSpec):
    def __init__(self):
        @composite
        def symbol_hyphen_weight(draw):
            sym = draw(sampled_from(el_symbols))
            weight = draw(integers(min_value=1, max_value=294))
            return f"{sym}-{weight}"

        @composite
        def weight_then_symbol(draw):
            sym = draw(sampled_from(el_symbols))
            weight = draw(integers(min_value=1, max_value=294))
            return f"{weight}{sym}"

        super().__init__(one_of(symbol_hyphen_weight(), weight_then_symbol())) # in principle, this'd be cleaner were it a class attribute.

    def phits(self, val):
        return val

    def python(self, val):
        return val

chemicals = {"H20", "CO2", "NH2", "NH3", "SF6", "TeF6", "CH4", "CH3", "C2H2", "C2H4", "C2H6", "C6H6", "CH32N3"}
class Chemical(ValSpec):
    def __init__(self):
        super().__init__(sampled_from(chemicals))

    def phits(self, val):
        if val in chemicals:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a valid chemical; got {val}."))

    def python(self, val):
        return val

class Orientation(ValSpec):
    def __init__(self):
        super().__init__(one_of(just("<"), just(">")))

    def phits(self, val):
        if val in ["<", ">"]:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be either `<` or `>`; got {val}."))

class Path(ValSpec):
    def __init__(self):
        super().__init__(text(whitelist_categories=("Nd", "L")))

    def phits(self, val):
        if isinstance(val, str):
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a filename string; got {val}."))

class Text(ValSpec):
    def __init__(self):
        super().__init__(text())

    def phits(self, val):
        if isinstance(val, str):
            return "{" + val + "}" # TODO: un-specialize if unnecessary
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a string; got {val}."))
