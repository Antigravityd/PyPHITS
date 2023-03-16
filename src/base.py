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

# Configuration options
g_value_type = "Python"
def settings(value_type="Python"):
    """Configure module-level parameters.

    Options:
      - value_type = "Python" | "PHITS"
        If "Python", use a remapped syntax that's more informative. E.g., `Parameters(control="output_echo_only")`
        instead of `Parameters(icntl=3)`.

        If "PHITS", use a syntax identical to that specified in the manual (i.e., the latter form in the above example).
        Note that some identifiers in PHITS do not conform to Python's identifier syntax (e.g. 2d-type, emin(14));
        these identifiers are sanitized as follows:
           - dashes -> underscores
           - beginning with a number -> that clause moved to the end
           - parentheses -> omitted

        For the examples above, the sanitized identifiers would be type_2d and emin14.
    """
    global g_value_type
    g_value_type = new_value_type


class PhitsObject:
    r"""The base class distinguishing objects that are intended to end up in some section of a .inp file,
    and defining equality and hashability of such objects in sensible ways.

    PhitsObject values correspond to some section of an input file.
    As such, they have a .definition() method that returns a textual representation of the value to be inserted into the input file:
    >>> print(Cylindrical("241Am", 2.2, fissile="neutrons", bounds=(-0.25, 0.25), r_out=0.3).definition())

    They also have a .prelude_str() method that returns text that is to preceed the whole set of values of that type:
    >>> print()

    Last, there's a .section_title() method that gives the name of the section into which the objects definition will be placed:
    >>> print()
n
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

    There's also a series of parameters related to grouping, that allow TODO: finish

    There's a "parser" and "validator" parameters, the first of which is a string representing part of a Lark parser
    used to extract the object from a string representing part of an .inp (it's not complete---for brevity,
    a set of common parse rules are retained by the this base class, and the parser string is appended to that)
    and the latter of which is a function which inspects the resulting tree and decides if it's legal PHITS input.
    """
    required = []
    positional = []
    optional = []
    ident_map = dict()
    value_map = dict()
    subobjects = []
    nones = dict()
    index = None
    no_hash = {"index", "value_map", "ident_map", "nones", "shape", "subobjects", "required", "positional", "optional",
               "group_by", "separator", "prelude", "max_groups", "group_size", "parser", "validator"}
    names = {"parameters", "source", "material", "surface", "cell", "transform", "temperature","mat_time_change","magnetic_field",
             "neutron_magnetic_field", "mapped_magnetic_field", "uniform_electromagnetic_field", "mapped_electromagnetic_field",
             "delta_ray", "track_structure", "super_mirror", "elastic_option", "importance", "weight_window", "ww_bias",
             "forced_collisions", "repeated_collisions", "volume", "multiplier", "mat_name_color", "reg_name", "counter", "timer",
             "t-track", "t-cross", "t-point", "t-adjoint", "t-deposit", "t-deposit2", "t-heat", "t-yield", "t-product", "t-dpa",
             "t-let", "t-sed", "t-time", "t-interact", "t-dchain", "t-wwg", "t-wwbg", "t-volume", "t-gshow", "t-rshow","t-3dshow"}

    def __init__(self, *args,  **kwargs):
        assert self.name in self.names, f"Unrecognized PHITS type {self.name} in PhitsObject initialization."

        if hasattr(self, syntax): # new-style definition, with high-resolution mapping
            required = list(map(lambda tup: tup[0],
                                sorted([(k, v) for k, v in self.syntax.items() if v[2] is not None], key=lambda tup: tup[1][2])))
            assert len(args) == len(required), f"Wrong number of positional arguments specified in the definition of {self.name} object."
            for idx, iden in enumerate(args):
                setattr(self, required[idx], arg if not isinstance(arg, list) else tuple(arg))

            for arg in self.syntax:
                if arg in kwargs:
                    setattr(self, arg, kwargs[arg] if not isinstance(arg, list) else tuple(arg))
                else:
                    setattr(self, arg, None)

            # TODO: reconsider
            for attr in self.subobjects:
                child = getattr(self, attr)
                if hasattr(child, self.name):
                    val = getattr(child, self.name)
                    if val is None:
                        setattr(child, self.name, self)

            remaining = {k: v for k, v in kwargs.items() if k not in self.syntax}
            if remaining:
                self.parameters = Parameters(**remaining)

        else: # Old-style definition, with "required" and "optional" lists, separate mappings, no value reassignment, etc.
              # Should be removed eventually.

            if len(args) == len(self.positional):
                for idx, arg in enumerate(args):
                    setattr(self, self.positional[idx], arg if not isinstance(arg, list) else tuple(arg))
            else:
                raise TypeError(f"Wrong number of positional arguments specified in the definition of {self.name} object.")

            for arg in self.required:
                if arg not in self.positional:
                    if arg in kwargs:
                        setattr(self, arg, kwargs[arg] if not isinstance(kwargs[arg], list) else tuple(kwargs[arg]))
                    else:
                        raise TypeError(f"Missing required argument in the definition of {self.name} object.")

            for arg in self.optional:
                if arg in kwargs:
                    setattr(self, arg, kwargs[arg] if not isinstance(kwargs[arg], list) else tuple(kwargs[arg]))
                else:
                    if arg in self.nones:
                        setattr(self, arg, self.nones[arg])
                    else:
                        setattr(self, arg, None)

            # TODO: reconsider
            for attr in self.subobjects:
                child = getattr(self, attr)
                if hasattr(child, self.name):
                    val = getattr(child, self.name)
                    if val is None:
                        setattr(child, self.name, self)



            remaining = {k: v for k, v in kwargs.items() if k not in self.required and k not in self.optional}
            if remaining:
                self.parameters = Parameters(**remaining)



    def add_definition(self, how, to, assignments=True):
        def idx(ob):
            if isinstance(ob, PhitsObject):
                return str(ob.index)
            else:
                return ob

        def attr_map(ob):
            if ob in self.ident_map:
                return self.ident_map[ob]
            else:
                return ob

        def val_map(ob):
            if ob in self.value_map:
                return self.value_map[ob]
            elif isinstance(ob, tuple):
                return " ".join(map(val_map, ob))
            elif isinstance(ob, Mesh): # TODO: consider refactoring out into a method of the Mesh class; it seems out of place here
                inp = f"{ob.axis[0]}-type = 1\n"
                inp += f"n{ob.axis[0]} = {len(ob.bins)-1}\n"
                for i in ob.bins:
                    inp += f"{i} "
                inp += "\n"
                return inp
            else:
                return ob
        
        for attr in how:
            if isinstance(attr, str) and attr[0] == "\\": # wat the fug
                to += attr[1:]
                return to
            assign = ""
            endstr = "\n" if assignments else " "
            spacing = " "

            if isinstance(attr, str) and attr[-1] == "\\":
                endstr = " "
                spacing = ""
                attr = attr[:-1]

            if callable(attr):
                to += attr(self)
                to += endstr

            elif isinstance(attr, tuple):
                to += self.add_definition(attr, to, assignments=False)
                if isinstance(attr[-1], str) and attr[-1][-1] == "\\":
                    to += " "
                else:
                    to += "\n"

            elif hasattr(self, attr):
                val = getattr(self, attr)
                if val is not None:
                    if attr in self.ident_map and isinstance(attr_map(attr), tuple):
                        for i, (att2, val2) in enumerate(zip(attr_map(attr), val)):
                            assign = f"{att2}{spacing}={spacing}" if assignments else ""
                            to += f"{assign}{idx(val_map(val2))}{endstr}"
                    else:

                        assign = f"{attr_map(attr)}{spacing}={spacing}" if assignments else ""
                        to += f"{assign}{idx(val_map(val))}{endstr}"

            else:
                if attr == "self":
                    to += idx(self)
                    to += endstr
                else:
                    to += attr
                    to += endstr
        return to

    def prelude_str(self):
        inp = self.add_definition(self.prelude, "")

        return continue_lines(inp)

    def definition(self):
        inp = self.add_definition(self.shape, "")

        return continue_lines(inp)

    def section_title(self):
        sec_name = self.name.replace("_", " ").title()
        return f"[{sec_name}]\n"




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
    def __init__(self, **kwargs):
        self.name = "parameters"
        for k, v in kwargs.items():
            assert k in self.syntax, f"Unrecognized parameter {k} (with value {v}) in initialization of Parameters object.
            Check that the correct parameters were passed to other objects."
            setattr(self, k, v)

    def __getitem__(self, key):
        return self.__dict__[key]

    def empty(self):
        return True if self.__dict__ == {"name": "parameters"} else False

    def definition(self):
        inp = ""
        for var, val in  self.__dict__.items():
            if var != "name":
                inp += f"{var} = {val}\n"

        return inp

