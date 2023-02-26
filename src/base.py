from collections import namedtuple, Iterable
from re import IGNORECASE
from lark import Lark, Transformer, Tree, Token
from lark.reconstruct import Reconstructor
from math import exp, log, log10, sqrt, acos, asin, atan, atan2, cos, cosh, sin, sinh, tan, tanh, pi
def continue_lines(inp):        # Continue lines that are too long
    r = ""
    length = 0
    last_whitespace = 0
    for i, char in enumerate(inp):
        if char.isspace():
            last_whitespace = i
        if length < 195:
            if char == "\n":
                length = 0
            r += char
            length += 1
        else:
            r = r[:last_whitespace] + " \\\n     " + r[last_whitespace:i] + char + r[i:]
            length = 0

    return r



class PhitsObject(Transformer):
    r"""The base class distinguishing objects that are intended to end up in some section of a .inp file,
    and defining equality and hashability of such objects in sensible ways.

    PhitsObject values correspond to some section of an input file.
    As such, they have a .definition() method that returns a textual representation of the value to be inserted into the input file:
    >>> print( Cylindrical("241Am", 2.2, fissile="neutrons", bounds=(-0.25, 0.25), r_out=0.3).definition() )

    They also have a .prelude() method that returns the line that is to preceed the whole set of values of that type:
    >>> print()

    Last, there's a .section_title() method that gives the name of the section into which the objects definition will be placed:
    >>> print()

    The PhitsObject class is also a factory for subtypes.
    As an example of how subtypes should be defined:
    >>> class Cylindrical(PhitsObject):
    ...    name = "source"
    ...    required = ["projectile", "energy"]
    ...    positional = ["projectile", "energy"]
    ...    optional = ["spin", "mask", "transform", "weight", "factor", "charge_override", "fissile", \
    ...                "center", "bounds", "r_out", "r_in", "elevation", "azimuth", "dispersion"]
    ...    ident_map = {"spin": ("sx", "sy", "sz"), "mask": ("reg", "ntmax"), "transform": "trcl", \
    ...                 "weight": "wgt", "charge_override": "izst", "fissile": "ispfs", "center": ("x0", "y0"), \
    ...                 "bounds": ("z0", "z1"), "r_out": "r0", "r_in": "r1", "elevation": "dir", "azimuth": "phi", \
    ...                 "dispersion": "dom", "energy": "e0", "projectile": "proj"}
    ...    value_map = {"neutrons": 2, True: 1}
    ...    shape = ("s-type = 1", "projectile", "spin", "mask", "transform", "weight", "factor", "charge_override", \
    ...             "fissile", "center", "bounds", "r_out", "r_in", "elevation", "azimuth", "dispersion", "energy")

    The "name" attribute is the [Section] of the input file into which the definition is to be inserted.

    "Required" gives the required (keyword or positional) arguments to the constructor, "positional" those which must be positional,
    and "optional" the optional (keyword or positional) arguments.

    The "ident_map" is a dictionary whose keys are arguments to the constructor (which become attributes of the instance),
    and whose values are the identifier(s) to which those arguments are to be converted for insertion into the input file.
    The inserted line is always something like "r0 = 5.3", hence the term "identifier" (the text to the left of the equality is substituted)
    This allows more idiomatic and descriptive naming of parameters, for example, "charge_override" as opposed to "izst" above,
    or the passage of iterables for assignments that ought to be grouped, like the components of a vector quantity.
    >>> print() # charge_override -> izst example
    >>> print() # spin -> sx, sy, sz example
    Similarly, the "value_map" does the same for the other side of the equals sign.

    Not shown in the above example, the "subobjects" parameter is used to indicate the names of any sections whose types can appear as attributes.
    These subobjects usually need a reference to the current object in their definition, so this parameter indicates to later processing
    to go in and update the subobjects accordingly
    This is used, for example, in the Material() object to enable the passage of TimeChange() objects directly to the material that is to change:
    >>> print()
    The "nones" parameter, also not shown above, is used to set a default value to a parameter that we consider optional, but that PHITS doesn't.
    >>> print()

    The real magic is in the "shape" parameter.
    It's a purpose-specific, ugly, and questionably-implemented analog of Emacs Lisp's skeleton system for programmatic text insertion,
    encoding how the attributes of the object are to be translated into text in the input file.
    The lion's share of data in input files are either parameter lines, of the form "<identifier> = <value>", or grid-like lines,
    of the form "<datum> <datum> <datum> ..."
    The value of the parameter is a tuple whose entries represent lines of input, ordered as given.
    If an entry is the name of an attribute of the object, name and value of that attribute are inserted as the parameter line "name = value".
    If an entry is a string that is not the name of an attribute, then it is inserted verbatim as text.
    If an entry is callable, it is considered to be a function that takes the current object as an argument and does something more complex.
    If an entry is a tuple itself, then this indicates a grid-like line, where the entries in the tuple are evaluated exactly as above,
    only if the entry is an attribute then only its value is inserted, and with spaces rather than newlines separating representations of
    entries of the tuple.
    Trailing backslashes in a string disable the insertion of the separating whitespace,
    and a string entry with only a backslash removes the separating whitespace of the previous entry (in case it's callable).

    There's also a series of parameters related to grouping, that allow objects to correspond to, say, rows in a grid in an input file.
    group_by gives the key on which to split these tables, prelude gives a setting (with syntax identical to shape) through which
    settings before the prelude can be altered, separator gives a string to be inserted between the groups (like a duplicate section header),
    max_groups indicates the limit on the number of groups (if one is set, say, by PHITS), and if group_size is set at construction,
    it will be reset to the size of a group an object is contained in when writing to input files.

    There're "parser" and "transformer" attributes, the first of which is a string representing part of a Lark parser
    used to extract the object from a string representing part of an .inp (it's not complete---for brevity,
    a set of common parse rules are retained by the this base class, and the parser string is appended to that)
    and the latter of which is a Lark transformer which turns the resulting parse tree into a value of the parent PhitsObject
    """
    name = "base"
    required = []
    positional = []
    optional = []
    ident_map = dict()
    value_map = dict()
    subobjects = []
    nones = dict()
    index = None
    no_hash = {"index", "value_map", "ident_map", "nones", "shape", "subobjects", "required", "positional", "optional",
               "group_by", "separator", "prelude", "max_groups", "group_size", "parser", "transformer"}
    names = {"base", "parameters", "source", "material", "surface", "cell", "transform", "temperature","mat_time_change","magnetic_field",
             "neutron_magnetic_field", "mapped_magnetic_field", "uniform_electromagnetic_field", "mapped_electromagnetic_field",
             "delta_ray", "track_structure", "super_mirror", "elastic_option", "importance", "weight_window", "ww_bias",
             "forced_collisions", "repeated_collisions", "volume", "multiplier", "mat_name_color", "reg_name", "counter", "timer",
             "t-track", "t-cross", "t-point", "t-adjoint", "t-deposit", "t-deposit2", "t-heat", "t-yield", "t-product", "t-dpa",
             "t-let", "t-sed", "t-time", "t-interact", "t-dchain", "t-wwg", "t-wwbg", "t-volume", "t-gshow", "t-rshow","t-3dshow"}
    def __init__(self, *args,  **kwargs):
        assert self.name in self.names, f"Unrecognized PHITS type {self.name} in PhitsObject initialization."

        if len(kwargs) == 0 and len(args) == 1 and args[0].isinstance(str): # initialization from a section of a PHITS .inp string
            return from_inp(args[0])
        else:
            for idx, arg in enumerate(self.required):
                if arg in kwargs:
                    setattr(self, arg, kwargs[arg] if not isinstance(kwargs[arg], list) else tuple(kwargs[arg]))
                elif len(args) >= idx + 1:
                    setattr(self, arg, args[idx] if not isinstance(args[idx], list) else tuple(args[idx]))
                else:
                    raise TypeError(f"Missing required argument {arg} in the definition of {self.name} object.")

            remaining = dict()
            for arg, val in kwargs.items():
                if arg in mapping:
                    setattr(self, arg, val)
                else:
                    remaining[arg] = val

            # TODO: rethink this
            for attr in self.subobjects:
                child = getattr(self, attr)
                if hasattr(child, self.name):
                    val = getattr(child, self.name)
                    if val is None:
                        setattr(child, self.name, self)

            if remaining:
                self.parameters = Parameters(**remaining)



    # this method is a candidate for a move
    def idx(ob):
        if isinstance(ob, PhitsObject):
            return str(ob.index)
        else:
            return ob

    def ident_to_phits(self, attr):
        if attr in self.mapping:
            return self.mapping[attr][0]
        else:
            return attr

    def ident_to_py(self, attr):
        rev = {v[0]: k for k, v in self.mapping.items()}
        if attr in rev:
            return rev[attr]
        else:
            return attr


    def value_to_py(self, attr, val):
        rev = {v[0]: k for k, v in self.mapping.items()}
        if attr in rev:
            res = self.mapping[rev[attr]][2](val, "Python")
            if callable(res):
                raise res(attr)
            else:
                return res
        else:
            return val

    def value_to_phits(self, attr, val):
        if attr in self.mapping:
            res = self.mapping[attr][2](val, "PHITS")
            if callable(res):
                raise res(attr)
            else:
                return res
        else:
            return val


    def definition(self):
        pars = Lark(self.parser + PhitsObject.g_parser, g_regex_flags=IGNORECASE, maybe_placeholders=False, ambiguity="resolve")
        return self.section_title() + continue_lines(Reconstructor(pars).reconstruct(self.tree()))

    def section_title(self):
        sec_name = self.name.replace("_", " ").title()
        return f"[{sec_name}]\n"

    g_parser = r"""
        IDENTIFIER.2: /[a-z0-9-]+/
        POSINT.5: /[0-9]+/
        INT.4: /-?[0-9]+/
        NUMBER.3: INT ["." [POSINT]] ["e" INT]
        PATHNAME: "/"? ((ARRAYENTRY | /[a-z][a-z0-9_.-]+/) "/")* (ARRAYENTRY | /[a-z][a-z0-9_.-]+/) "/"?
        ARRAYENTRY: IDENTIFIER "(" POSINT ")"
        function{name}: name "(" computation ")"
        binop{sym}: computation sym computation
        PI.1: /pi/
        VARIABLE: /x/
        computation: " "* (NUMBER | PI
                            | binop{/\+/} | binop{/-/} | binop{/\*/} | binop{/\//} | binop{/\*\*/}
                            | function{/float/} | function{/int/} | function{/abs/} | function{/exp/} | function{/log/} | function{/log10/}
	      	            | function{/max/} | function{/min/} | function{/mod/} | function{/nint/} | function{/sign/} | function{/sqrt/}
	      	            | function{/acos/} | function{/asin/} | function{/atan/} | function{/atan2/} | function{/cos/}
	      	            | function{/cosh/} | function{/sin/} | function{/sinh/} | function{/tan/} | function{/tanh/}) " "*
        functiondef: " "* (NUMBER | PI | VARIABLE
                           | binop{/\+/} | binop{/-/} | binop{/\*/} | binop{/\//} | binop{/\*\*/}
                           | function{/float/} | function{/int/} | function{/abs/} | function{/exp/} | function{/log/} | function{/log10/}
	      	           | function{/max/} | function{/min/} | function{/mod/} | function{/nint/} | function{/sign/} | function{/sqrt/}
	      	           | function{/acos/} | function{/asin/} | function{/atan/} | function{/atan2/} | function{/cos/}
	      	           | function{/cosh/} | function{/sin/} | function{/sinh/} | function{/tan/} | function{/tanh/}) " "*
        escapedcomputation: "{" computation "}" -> computation
        numbergrid: numberline+
        numberline: " "* (computation " "+)+ computation? "\n"
        columntitle: " "* (IDENTIFIER " "+)+ IDENTIFIER? "\n"
        %ignore /\n(?=\n)/
        %import common.WS
        """
    IDENTIFIER = str
    POSINT = int
    INT = int
    NUMBER = float
    FILENAME = str
    PATHNAME = str
    ARRAYENTRY = str
    PI = pi
    def computation(self, comp):

        if not isinstance(comp[0], Tree):
            return comp[0]
        elif comp[0].data == "function":
            matched = comp[0].children
            fname = matched[0].lower()
            try:
                return eval(fname)(matched[1])
            except NameError: # function not immediately Pythonizable
                if fname == "nint":
                    return round(matched[1])
        elif comp[0].data == "binop":
            matched = comp[0].children
            return eval(f'lambda x, y: x {matched[1]} y')(matched[0], matched[2])

    numbergrid = list
    numberline = list
    # TODO: consider whether functiondef deserves conversion into Python function object or not; we'd have to convert back which may suck.


    @classmethod
    def from_inp(self, segment):
        parser = Lark(self.parser + PhitsObject.g_parser, g_regex_flags=IGNORECASE)
        transformer = PhitsObject() * self()
        tree = parser.parse(segment)
        breakpoint()
        return transformer.transform(tree)


    def __eq__(self, other):
        if type(self) != type(other):
            return False
        elif hasattr(self, "__dict__") and hasattr(other, "__dict__"):
            d1 = {k: v for k, v in self.__dict__.items() if k not in self.no_hash}
            d2 = {k: v for k, v in other.__dict__.items() if k not in self.no_hash}
            return d1 == d2

    def __hash__(self):
        return hash(tuple(v for k, v in sorted(self.__dict__.items()) \
                          if (self not in v.__dict__.values() if hasattr(v, "__dict__") else True) and k not in self.no_hash))



# Some helper functions useful in consisely defining value_map parameters.
# The non-inv functions go from Python -> PHITS format, and the inv's go the opposite.
# Each function must return a function that validates the input value of a given parameter and outputs a valid version of it in the output;
# e.g. the non-inv functions are expected to be called on values of parameters passed via Python, validate they're acceptable,
# and produce the corresponding parameter value that could be pasted into the PHITS file.
# If a parameter fails validation, a descriptive ValueError should be returned.
# The PHITS -> Python is guaranteed to have the correct type, conditional on parser correctness, but Python input needs checked.


def choice10(c_style=False, true=True, false=False):
    def r(val, out):
        if out == "PHITS":
            if val == true:
                return 0 if c_style else 1
            elif val == false:
                return 1 if c_style else 0
            else:
                return lambda var: ValueError(f"`{var}` must be either True or False; got {val}.")
        elif out == "Python":
            if val == 0:
                return true if c_style else false
            elif val == 1:
                return false if c_style else true
            else:
                return lambda var: ValueError(f"`{self.ident_map[var]}` must be either 0 or 1; got {val}.")
    return r



def integer(val, out):
    if out == "PHITS":
        if isinstance(val, int):
            return val
        else:
            return lambda var: ValueError(f"`{var}` must be an integer; got {val}.")
    elif out == "Python":
        if x % 1 == 0:
            return val
        else:
            return lambda var: ValueError(f"`{self.ident_map[var]}` must be an integer; got {val}.")



def real(val, out):
    if out == "PHITS":
        if isinstance(val, float):
            return val
        else:
            return lambda var: ValueError(f"`{var}` must be a float; got {val}.")
    elif out == "Python":
        return val



def posint(val, out):
    if out == "PHITS":
        if isinstance(val, int) and val > 0:
            return val
        else:
            return lambda var: ValueError(f"`{var}` must be a positive integer; got {val}.")
    elif out == "Python":
        if val > 0:
            return val
        else:
            return lambda var: ValueError(f"`{self.ident_map[var]}` must be positive; got {val}.")



def posreal(val, out):
    if out == "PHITS":
        if isinstance(val, float) and val > 0:
            return val
        else:
            lambda var: ValueError(f"`{var}` must be a positive floating-point value; got {val}.")
    elif out == "Python":
        if val > 0:
            return val
        else:
            lambda var: ValueError(f"`{self.ident_map[var]}` must be positive; got {val}.")





def triplechoice(val, out):
    if out == "PHITS":
        if isinstance(val, int):
            return True if val > 0 else False
        else:
            return lambda var: ValueError(f"`{var}` must be an integer; got {val}.")
    elif out == "Python":
        if val > 0:
            return True
        else:
            return False








def negdisable(val, out):
    if out == "PHITS":
        if isinstance(val, int) and val > 0:
            return val
        elif val is None:
            return -1.0
        else:
            return lambda var: ValueError(f"`{self.ident_map[var]}` must be a positive integer or None; got {val}.")
    elif out == "Python":
        if val > 0:
            return val
        else:
            return None




def between(start, stop):
    def r(val, out):
        if out == "PHITS":
            if isinstance(val, int) and val >= start and val <= stop:
                return val
            else:
                return lambda var: ValueError(f"`{var}` must be an integer between {start} and {stop}, inclusive; got {val}.")
        elif out == "Python":
            if val >= start and val <= stop:
                return val
            else:
                return lambda var: ValueError(f"`{self.ident_map[var]}` must be between {start} and {stop}, inclusive; got {val}.")
    return r



def zerospecial(val, out):
    if out == "PHITS":
        if isinstance(val, int):
            if val == 0:
                return zero
            else:
                return val
        else:
            return lambda var: ValueError(f"`{var}` must be an integer; got {val}.")
    elif out == "Python":
        if val == zero:
            return 0
        else:
            return val
    return r

# def allornothing():
#     def r(val, out):
#         if out == "PHITS":
#             if

def finbij(dic):
    def r(val, out):
        if out == "PHITS":
            if val in dic:
                return dic[val]
            else:
                return lambda var: ValueError(f"`{var}` must be one of {list(dic.keys())}; got {val}.")
        elif out == "Python":
            rev = {v: k for k, v in dic.items()}
            if val in rev:
                return rev[val]
            else:
                return lambda var: ValueError(f"`{self.ident_map[var]}` must be one of {list(rev.keys())}; got {val}.")


    return r

class Parameters(PhitsObject, Transformer):
    """A "dictionary with an attitude" representing an entry in the [Parameters] section of an input file.
    Any extra keyword arguments to any constructors are minted into parameter objects.


    >>> print(Parameters(ndedx=2, dbcutoff=3.3).definition())
    ndedx = 2
    dbcutoff = 3.3

    """
    name = "parameters"
    # implicitly defines the __init__ and what arguments are acceptable to it
    # "Python argument name": ("PHITS name", default, fn<validate_and_remap>)
    mapping = {"control": ("icntl", "normal",
                           finbij({"normal": 0, "output_cross-section": 1, "output_echo_only": 3, "all_reg_void": 5,
                                   "source_check": 6, "show_geometry": 7, "show_geometry_with_xyz": 8, "show_regions": 9,
                                   "show_regions_with_tally": 10, "show_3d_geometry": 11, "use_dumpall": 12, "sum_tally": 13,
                                   "auto_volume": 14, "ww_bias_tally": 15, "analysis_script": 16, "anatally": 17})),
               "max_histories": ("maxcas", 10, posint),
               "max_batches": ("maxbch", 10, posint),
               "nuclear_memory_rescale": ("xsmemory", 1.0, posreal),
               "timeout": ("timeout", None, negdisable),
               "stdev_control": ("istdev", "normal", finbij({-2:"history_restart", -1:"batch_restart", 0:"normal", 1:"batch", 2:"history"})),
               "share_tallies": ("italsh", False, choice10()),
               "check_consistency": ("ireschk", False, choice10(c_style=True)),
               "xor_prng": ("nrandgen", True, choice10()),
               "seed_skip": ("irskeep", 0, integer),
               "random_seed": ("rseed", 6.657299061401e12, real),
               "seed_from_time": ("itimrand", False, choice10()),
               # bitrseed???
               "proton_e_cutoff": ("emin(1)", 1.0e-3, posreal),
               "neutron_e_cutoff": ("emin(2)", 1.0e-11, posreal),
               "pion+_e_cutoff": ("emin(3)", 1.0e-3, posreal),
               "pion0_e_cutoff": ("emin(4)", 1.0e-3, posreal),
               "pion-_e_cutoff": ("emin(5)", 1.0, posreal),
               "muon+_e_cutoff": ("emin(6)", 1.0e-3, posreal),
               "muon-_e_cutoff": ("emin(7)", 1.0e-3, posreal),
               "kaon+_e_cutoff": ("emin(8)", 1.0e-3, posreal),
               "kaon0_e_cutoff": ("emin(9)", 1.0e-3, posreal),
               "kaon-_e_cutoff": ("emin(10)", 1.0e-3, posreal),
               "other_e_cutoff": ("emin(11)", 1.0, posreal),
               "electron_e_cutoff": ("emin(12)", 1.0e9, posreal),
               "positron_e_cutoff": ("emin(13)", 1.0e9, posreal),
               "photon_e_cutoff": ("emin(14)", 1.0e-3, posreal),
               "deuteron_e_cutoff": ("emin(15)", 1.0e-3, posreal),
               "triton_e_cutoff": ("emin(16)", 1.0e-3, posreal),
               "he-3_e_cutoff": ("emin(17)", 1.0e-3, posreal),
               "he-4_e_cutoff": ("emin(18)", 1.0e-3, posreal),
               "nucleon_e_cutoff": ("emin(19)", 1.0e-3, posreal),
               "proton_e_max": ("dmax(1)", 1.0e-3, posreal),
               "neutron_e_max": ("dmax(2)", 20.0, posreal),
               "electron_e_max": ("dmax(12)", 1.0e9, posreal),
               "positron_e_max": ("dmax(13)", 1.0e9, posreal),
               "photon_e_max": ("dmax(14)", 1.0e3, posreal),
               "deuteron_e_max": ("dmax(15)", 1.0e-3, posreal),
               "he-4_e_max": ("dmax(18)", 1.0e-3, posreal),
               "photonuclear_e_max": ("dpnmax", 0.0, posreal),
               # lib(i)
                "charged_e_min": ("esmin", 0.001, posreal),
               "charged_e_max": ("esmax", 300000.0, posreal),
               # cmin
                "electron-positron_track-structure_e_min": ("etsmin", 1e-6, posreal),
               "electron-positron_track-structure_e_max": ("etsmax", 1e-2, posreal),
               "nucleon_track-structure_e_max": ("tsmax", 1e-3, posreal),
               "electric_transport_type": ("negs", "PHITS", finbij({"PHITS": -1, "ignore": 0, "EGS5": 1})),
               "automatic_e_bounds": ("nucdata", True, choice10()),
               "electron-positron_adjust_weight_over_e_max": ("ieleh", False, choice10()),
               "nucleon-nucleus_model_switch_e": ("ejamnu", 20.0, posreal),
               "pion-nucleus_model_switch_e": ("ejampi", 20.0, posreal),
               "isobar_max_e": ("eisobar", 0.0, posreal),
               "isobar_model": ("isobar", False, choice10()),
               # etc.
                }
    implications = "" # the value of some options may restrict the value of others. ("ident1", fn<pred>): ("ident2": fn<pred>)


    # parser/transformer stuff
    parser = r"""
        %ignore WS
        start: assignment*
        assignment: iden "=" val
        iden: IDENTIFIER | ARRAYENTRY
        val: computation | PATHNAME"""

    def iden(self, ide):
        return str(ide[0])

    def val(self, va):
        return float(va[0])

    def assignment(self, ass):
        return (self.ident_to_py(ass[0]), self.value_to_py(ass[0], ass[1]))

    def start(self, asses):
        return Parameters(**dict(asses))

    def tree(self):
        start = []
        for k, v in self.__dict__.items(): # break out into cases to handle more complicated remapping
            outi = self.ident_to_phits(k)
            outv = self.value_to_phits(k, v)
            if "(" in outi:
                if isinstance(outv, str):
                    start.append(Tree(Token('RULE', 'assignment'),
                                      [Tree(Token("RULE", 'iden'), [Token("ARRAYENTRY", outi)]),
                                       Tree(Token("RULE", 'val'), [Token("PATHNAME", outv)])]))
                elif isinstance(outv, float):
                    start.append(Tree(Token('RULE', 'assignment'),
                                      [Tree(Token("RULE", 'iden'), [Token("ARRAYENTRY", outi)]),
                                       Tree(Token("RULE", 'val'), [Tree(Token("RULE", "computation"), [outv])])]))
            else:
                if isinstance(outv, str):
                    start.append(Tree(Token('RULE', 'assignment'),
                                      [Tree(Token("RULE", "iden"), [Token("IDENTIFIER", outi)]),
                                       Tree(Token("RULE", 'val'), [Token("PATHNAME", outv)])]))
                elif isinstance(outv, float) or isinstance(outv, int):
                    start.append(Tree(Token('RULE', 'assignment'),
                                      [Tree(Token("RULE", 'iden'), [Token("IDENTIFIER", outi)]),
                                       Tree(Token("RULE", 'val'), [Tree(Token("RULE", "computation"), [outv])])]))
        r =  Tree(Token('RULE', 'start'), start)
        breakpoint()

        return r


    # bookkeeping
    def __init__(self, **kwargs):
        for attr, val in kwargs.items():
            if attr in self.mapping:
                setattr(self, attr, val)
            else:
                raise ValueError( \
f"Argument `{attr}` (value `{val}`) not a valid global PHITS parameter---check that all object initializer keword arguments are correct.")


    def __getitem__(self, key):
        return self.__dict__[key]
    def empty(self):
        return True if self.__dict__ == {"name": "parameters"} else False


# I want to get rid of this
# class Mesh():
#     """Represents all list-typed data in PHITS."""
#     def __init__(self, axis, bins=None): # All bin generation is easily done in Python via list comprehensions
#         assert axis in ["energy", "time", "x", "y", "z", "radius", "angle", "let"], f"Unrecognized axis {axis} in mesh definition."
#         self.axis = axis
#         self.bins = tuple(bins)
#         print(self.bins)
#         if axis != "angle":
#             self.type = 2
#         else:
#             pass  # TODO: figure out angle mesh


#     def __eq__(self, other):
#         if type(self) != type(other):
#             return False

#         else:
#             return {k: v for k, v in self.__dict__.items() if v is not other} \
#                 == {k: v for k, v in other.__dict__.items() if v is not self}

#     def __hash__(self):
#         return hash(tuple(v for k, v in sorted(self.__dict__.items()) \
#                           if (self not in v.__dict__.values() if hasattr(v, "__dict__") else True)))
