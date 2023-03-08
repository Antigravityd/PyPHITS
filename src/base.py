from collections import Iterable
import regex as re
import dataclasses
from lark import Lark, Transformer, Tree, Token
from lark.reconstruct import Reconstructor
from math import exp, log, log10, sqrt, acos, asin, atan, atan2, cos, cosh, sin, sinh, tan, tanh, pi
from valspec import *
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

    There're "grammar" and "transformer" attributes, the first of which is a string representing part of a Lark grammar
    used to extract the object from a string representing part of an .inp (it's not complete---for brevity,
    a set of common parse rules are retained by the this base class, and the grammar string is appended to that)
    and the latter of which is a Lark transformer which turns the resulting parse tree into a value of the parent PhitsObject
    """
    name = "base"
    subobjects = []
    index = None
    no_hash = {"index", "value_map", "ident_map", "nones", "shape", "subobjects", "required", "positional", "optional",
               "group_by", "separator", "prelude", "max_groups", "group_size", "grammar", "transformer"}
    names = {"base", "parameters", "source", "material", "surface", "cell", "transform", "temperature","mat_time_change","magnetic_field",
             "neutron_magnetic_field", "mapped_magnetic_field", "uniform_electromagnetic_field", "mapped_electromagnetic_field",
             "delta_ray", "track_structure", "super_mirror", "elastic_option", "importance", "weight_window", "ww_bias",
             "forced_collisions", "repeated_collisions", "volume", "multiplier", "mat_name_color", "reg_name", "counter", "timer",
             "t-track", "t-cross", "t-point", "t-adjoint", "t-deposit", "t-deposit2", "t-heat", "t-yield", "t-product", "t-dpa",
             "t-let", "t-sed", "t-time", "t-interact", "t-dchain", "t-wwg", "t-wwbg", "t-volume", "t-gshow", "t-rshow","t-3dshow"}
    def __init__(self, *args,  **kwargs):
        # assert self.name in self.names, f"Unrecognized PHITS type {self.name} in PhitsObject initialization."

        if len(args) == 0 and len(kwargs) == 0: # empty initialization, probably for transformer instance
            return
        elif len(kwargs) == 0 and len(args) == 1 and isinstance(args[0], str): # initialization from a section of a PHITS .inp string
            for k, v in self.from_inp(args[0]).__dict__.items():
                setattr(self, k, v)
        else:
            required = [i for i in self.syntax if self.syntax[i][1] is None]   # TODO: think harder about this
            for idx, arg in enumerate(required):
                if arg in kwargs:
                    setattr(self, arg, kwargs[arg] if not isinstance(kwargs[arg], list) else tuple(kwargs[arg]))
                elif len(args) >= idx + 1:
                    setattr(self, arg, args[idx] if not isinstance(args[idx], list) else tuple(args[idx]))
                else:
                    raise TypeError(f"Missing required argument {arg} in the definition of {self.name} object.")

            remaining = dict()
            for arg, val in kwargs.items():
                if arg in self.syntax:
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

    @classmethod
    def inv_syntax(self):
        r = dict()
        for py_ident, (phits_ident, valspec, argorder) in self.syntax.items():
            if isinstance(phits_ident, tuple):
                for i, phits_ident_real in enumerate(phits_ident):
                    r[phits_ident_real] = (py_ident, valspec[i], argorder)
            else:
                r[phits_ident] = (py_ident, valspec, argorder)

        return r

    @staticmethod
    def sanitize(phits_iden):
        """Takes un-Larkizable names like '2d-type' and 'emin(2)' to alternatives.
        Does not need to be bidirectional, just consistent; both the grammar and the corresponding transformer functions are constructed,
        so this process merely need be univalent between those two constructions."""
        san = phits_iden.replace("(", "")
        san = san.replace(")", "")
        san = san.replace("-", "_")
        return san

    @classmethod
    def grammar(self):
        """Constructs and returns a full (save for global definitions in PhitsObject.g_grammar) grammar from the syntax self.syntax
        and the macro-ified grammar self._grammar.

        Macros:
          @assign_among|<iterable>| -- when <iterable> is a subset of the inv_syntax() dictionary, match assignment statements in that subset.
          @parse_of|<PhitsObject>|  -- match any valid definition of the PhitsObject in question, according to its grammar."""
        grammar = self._grammar
        def grammarmod(directive, replace, suffix=None):
            """Takes the grammar, finds all appearences of 'directive', calculates and performs textual replacement of the result of 'replace'
            (which is a fn syntax the result of the directive and match # to the substitution), and optionally appends the result of 'suffix'
            (a lambda syntax the list of all arguments found to all directives of the given type to a string) to the end of the grammar."""
            nonlocal grammar, self
            matches = re.finditer(f"(?<=@{directive})\\|(?:[^|]+|(?R))*+\\|", grammar)
            rules = []
            sub = []
            for i, match in enumerate(matches):
                scope = locals()
                arg = eval(match[0][1:-1], scope)
                sub = replace(i, arg)
                grammar = re.sub(f"@{directive}\\|(?:[^|]+|(?R))*+\\|", sub, grammar, count=1)
                if suffix is not None:
                    grammar += suffix(i, arg)



        # @assign_among
        def ass_repl(i, keys):
            nonlocal self
            inv = self.inv_syntax()
            subset = {k: inv[k] for k in keys}
            sub = []
            for phits_iden, (py_iden, val_spec, arg_order) in subset.items():
                san = self.sanitize(phits_iden)
                sub.append(san)

            return "(" + ' | '.join(sub) + ")"

        def ass_suff(i, subset):
            nonlocal self
            rules = []
            for phits_iden, (py_iden, val_spec, argorder) in subset.items():
                san = self.sanitize(phits_iden)
                other_idens = self.syntax[py_iden][0]
                rules.append(f"{san}: \"{phits_iden}\" \"=\" {val_spec.parse_rule} \"\\n\"")
            return "\n".join(rules)

        grammarmod("assign_among", ass_repl, ass_suff)

        # @parse_of
        def parse_repl(i, po): # NOTE: must be careful to avoid name collisions in subrules.
                               # It's possible to uniquify, of course, but too much work.
            nonlocal self
            return self.sanitize(po.name) + "start"


        def parse_suff(i, po):
            nonlocal self
            return re.sub("start:", self.sanitize(po.name) + "start:", po._grammar)


        grammarmod("parse_of", parse_repl, parse_suff)

        # @assign_then_grid
        def grid_repl(i, entr):
            nonlocal self
            sub = f"{self.sanitize(entr)}"
            return "(" + sub + ")"

        def grid_suff(i, entr):
            nonlocal self
            phits_iden = entr
            val_spec = self.inv_syntax()[entr][1]
            rule = f"""
            {self.sanitize(phits_iden)}: "{phits_iden}" "=" {val_spec.parse_rule} "\\n" numbergrid?
            """
            return rule

        grammarmod("assign_then_grid", grid_repl, grid_suff)


        return grammar + PhitsObject.g_grammar

    def definition(self):
        pars = Lark(self.grammar(), g_regex_flags=re.IGNORECASE, maybe_placeholders=False, ambiguity="resolve")
        return continue_lines(Reconstructor(pars).reconstruct(self.tree()))

    def section_title(self):
        sec_name = self.name.replace("_", " ").title()
        return f"[{sec_name}]\n"

    g_grammar = r"""
        IDENTIFIER.2: /[a-z0-9-]+/
        POSINT.5: /[1-9][0-9]*/
        INT.4: /-?/ (/0/ | POSINT)
        BIN: /0/ | /1/
        NUMBER.3: INT ["." [POSINT]] ["e" INT]
        PATHNAME: "/"? ((ARRAYENTRY | /[a-z][a-z0-9_.-]+/) "/")* (ARRAYENTRY | /[a-z][a-z0-9_.-]+/) "/"?
        ARRAYENTRY: IDENTIFIER "(" POSINT ")"
        function{name}: name "(" computation ")"
        binop{sym}: computation sym computation
        PI.1: /pi/
        VARIABLE: /x/
        computation: " "* (NUMBER | PI | VARIABLE
                            | binop{/\+/} | binop{/-/} | binop{/\*/} | binop{/\//} | binop{/\*\*/}
                            | function{/float/} | function{/int/} | function{/abs/} | function{/exp/} | function{/log/} | function{/log10/}
	      	            | function{/max/} | function{/min/} | function{/mod/} | function{/nint/} | function{/sign/} | function{/sqrt/}
	      	            | function{/acos/} | function{/asin/} | function{/atan/} | function{/atan2/} | function{/cos/}
	      	            | function{/cosh/} | function{/sin/} | function{/sinh/} | function{/tan/} | function{/tanh/}) " "*

        escapedcomputation: "{" computation "}" -> computation
        numbergrid: numberline+
        numberline: " "* (computation " "+)+ computation? "\n"
        columntitle: " "* (IDENTIFIER " "+)+ IDENTIFIER? "\n"
        SPACE: /[ \t]+/
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
        if isinstance(comp[0], Tree):
            if comp[0].data == "function":
                fname = comp[0].children[0]
                fname = "round" if fname == "nint" else fname
                arg = comp[0].children[1]
                if arg == "x":
                    return lambda x: eval(fname)(x)
                elif callable(arg):
                    return lambda x: eval(fname)(arg(x))
                else:
                    return eval(f"{fname}(arg)")
            elif comp[0].data == "binop":
                op = comp[0].children[1]
                arg1 = comp[0].children[0]
                arg2 = comp[0].children[2]
                unfix = eval(f"lambda x, y: x {op} y")
                if arg1 == "x" or arg2 == "x":
                    return eval(f"lambda x: {arg1} {op} {arg2}")
                elif callable(arg1) and not callable(arg2):
                    return lambda x: unfix(arg1(x), arg2)
                elif not callable(arg1) and callable(arg2):
                    return lambda x: unfix(arg1, arg2(x))
                elif callable(arg1) and callable(arg2):
                    return lambda x: unfix(arg1(x), arg2(x))
                else:
                    return unfix(arg1, arg2)
        else:
            return comp[0]

        
    numbergrid = list
    numberline = list


    @classmethod
    def from_inp(self, segment):
        transformer = Transformer()
        transformer.parsing = self
        for phits_ident, (py_ident, valspec, argorder) in self.inv_syntax().items():
            scope = locals() | {"syntax": self.syntax, "py_ident": py_ident, "phits_ident": phits_ident, "valspec": valspec,
                                "multiple": isinstance(self.syntax[py_ident][0], tuple)}

            # this eval hack is motivated by a local function definition simply...not working---the function got overwritten each iteration
            san = self.sanitize(phits_ident)
            exec(f"""def {san}(s, tr):
    if multiple:
        return (py_ident, valspec.python(tr[0]), syntax[py_ident][0].index(phits_ident))
    else:
        return (py_ident, valspec.python(tr[0]))""", scope)

            setattr(transformer, san, scope[san].__get__(transformer))


        def start(s, tr):
            if isinstance(tr[0], list):
                tr = tr[0]
            divided = dict()
            for tup in tr:
                py_ident = tup[0]
                py_val = tup[1]
                if len(tup) == 3:
                    argindex = tup[2]
                    if py_ident in divided:
                        divided[py_ident] |= {argindex: py_val}
                    else:
                        divided[py_ident] = {argindex: py_val}
                else:
                    divided[py_ident] = py_val

            for k, v in divided.items():
                if isinstance(v, dict):
                    divided[k] = tuple((v[i] for i in sorted(v.keys())))

            return s.parsing(**divided)

        setattr(transformer, "start", start.__get__(transformer))

        # Any methods in subclass not escaped with an underscore are added to transformer
        for attr in set(dir(self)) - set(dir(PhitsObject)):
            if callable(getattr(self, attr)) and attr[0] != "_":
                setattr(transformer, attr, getattr(self, attr).__get__(transformer))

        transformer = PhitsObject() * transformer

        grammar = Lark(self.grammar(), g_regex_flags=re.IGNORECASE)
        tree = grammar.parse(segment)
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







class Parameters(PhitsObject):
    """A "dictionary with an attitude" representing an entry in the [Parameters] section of an input file.
    Any extra keyword arguments to any constructors are minted into parameter objects.


    >>> print(Parameters(ndedx=2, dbcutoff=3.3).definition())
    ndedx = 2
    dbcutoff = 3.3
    """
    name = "parameters"
    # implicitly defines the __init__, what arguments are acceptable to it, and how to generate the parser and tree.
    # "Python argument name": ("PHITS name", default, subclass(ValSpec), RequirementGroup)
    syntax = {"control": ("icntl",
                           FinBij({"normal": 0, "output_cross-section": 1, "output_echo_only": 3, "all_reg_void": 5,
                                   "source_check": 6, "show_geometry": 7, "show_geometry_with_xyz": 8, "show_regions": 9,
                                   "show_regions_with_tally": 10, "show_3d_geometry": 11, "use_dumpall": 12, "sum_tally": 13,
                                   "auto_volume": 14, "ww_bias_tally": 15, "analysis_script": 16, "anatally": 17}), None),
               "max_histories": ("maxcas", PosInt(), None),
               "max_batches": ("maxbch", PosInt(), None),
               "nuclear_memory_rescale": ("xsmemory", 1.0, PosReal(), None),
               "timeout": ("timeout", NegDisable(), None),
               "stdev_control": ("istdev", FinBij({-2: "history_restart", -1: "batch_restart", 0: "normal", 1: "batch", 2: "history"}), None),
               "share_tallies": ("italsh", Choice10(), None),
               "check_consistency": ("ireschk", Choice10(c_style=True), None),
               "xor_prng": ("nrandgen", Choice10(), None),
               "seed_skip": ("irskeep", Integer(), None),
               "random_seed": ("rseed", Real(), None),
               "seed_from_time": ("itimrand", Choice10(), None),
               # bitrseed?,
               "proton_e_cutoff": ("emin(1)", PosReal(), None),
               "neutron_e_cutoff": ("emin(2)", PosReal(), None),
               "pionp_e_cutoff": ("emin(3)", PosReal(), None),
               "pion0_e_cutoff": ("emin(4)", PosReal(), None),
               "pionm_e_cutoff": ("emin(5)", PosReal(), None),
               "muonp_e_cutoff": ("emin(6)", PosReal(), None),
               "muonm_e_cutoff": ("emin(7)", PosReal(), None),
               "kaonp_e_cutoff": ("emin(8)", PosReal(), None),
               "kaon0_e_cutoff": ("emin(9)", PosReal(), None),
               "kaonm_e_cutoff": ("emin(10)", PosReal(), None),
               "other_e_cutoff": ("emin(11)", PosReal(), None),
               "electron_e_cutoff": ("emin(12)", PosReal(), None),
               "positron_e_cutoff": ("emin(13)", PosReal(), None),
               "photon_e_cutoff": ("emin(14)", PosReal(), None),
               "deuteron_e_cutoff": ("emin(15)", PosReal(), None),
               "triton_e_cutoff": ("emin(16)", PosReal(), None),
               "he3_e_cutoff": ("emin(17)", PosReal(), None),
               "he4_e_cutoff": ("emin(18)", PosReal(), None),
               "nucleon_e_cutoff": ("emin(19)", PosReal(), None),
               "proton_e_max": ("dmax(1)", PosReal(), None),
               "neutron_e_max": ("dmax(2)", PosReal(), None),
               "electron_e_max": ("dmax(12)", PosReal(), None),
               "positron_e_max": ("dmax(13)", PosReal(), None),
               "photon_e_max": ("dmax(14)", PosReal(), None),
               "deuteron_e_max": ("dmax(15)", PosReal(), None),
               "he4_e_max": ("dmax(18)", PosReal(), None),
               "photonuclear_e_max": ("dpnmax", PosReal(), None),
               # lib(i)
               "charged_e_min": ("esmin", PosReal(), None),
               "charged_e_max": ("esmax", PosReal(), None),
               # cmin
                "electron_positron_track_structure_e_min": ("etsmin", PosReal(), None),
               "electron_positron_track_structure_e_max": ("etsmax", PosReal(), None),
               "nucleon_track_structure_e_max": ("tsmax", PosReal(), None),
               "electric_transport_type": ("negs", FinBij({"PHITS": -1, "ignore": 0, "EGS5": 1}), None),
               "automatic_e_bounds": ("nucdata", Choice10(), None),
               "electron_positron_adjust_weight_over_e_max": ("ieleh", Choice10(), None),
               "nucleon_nucleus_model_switch_e": ("ejamnu", PosReal(), None),
               "pion_nucleus_model_switch_e": ("ejampi", PosReal(), None),
               "isobar_max_e": ("eisobar", PosReal(), None),
               "isobar_model": ("isobar", Choice10(), None),
               # etc.
                }
    implications = "" # the value of some options may restrict the value of others. ("ident1", fn<pred>): ("ident2": fn<pred>)

    def tree(self):
        tree = []
        for py_attr, py_val in self.__dict__.items():
            if py_attr in self.syntax:
                phits_attr = self.syntax[py_attr][0]
                valspec = self.syntax[py_attr][2]
                phits_val = valspec.phits(py_val)
                if "(" in phits_attr or ")" in phits_attr:
                    phits_attr = Token("ARRAYENTRY", phits_attr)
                else:
                    phits_attr = Token("IDENTIFIER", phits_attr)

                if isinstance(phits_val, str):
                    tree.append(Tree(Token("RULE", "assignment"), [phits_attr, Token("PATHNAME", py_val)]))
                else:
                    tree.append(Tree(Token("RULE", "assignment"), [phits_attr, Tree(Token("RULE", "computation"),
                                                                                    [Token("NUMBER", py_val)])]))

        return Tree(Token("RULE", "start"), tree)



    # grammar/transformer stuff
    _grammar = r"""
    %ignore SPACE
    start: assignment+
    assignment: @assign_among|self.inv_syntax().keys()|
    """

    def assignment(self, tr):
        return (self.inv_syntax()[tr[0]][0], self.inv_syntax()[tr[0]][2].python(tr[1]))

    def start(self, tr):
        return Parameters(**dict(tr))

    # bookkeeping
    def __init__(self, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], str):
            for k, v in self.from_inp(args[0]).__dict__.items():
                setattr(self, k, v)
        elif len(args) == 0:
            for attr, val in kwargs.items():
                if attr in self.syntax:
                    self.syntax[attr][2].phits(val)
                    setattr(self, attr, val)
                else:
                    raise ValueError( \
f"Argument `{attr}` (value `{val}`) not a valid global PHITS parameter---check that all object initializer keword arguments are correct.")
        else:
            raise ValueError(f"Parameters() takes exclusively keyword arguments unless initializing from .inp string; got {args}.")


    def __getitem__(self, key):
        return self.__dict__[key]

    def empty(self):
        return True if self.__dict__.keys() else False


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
