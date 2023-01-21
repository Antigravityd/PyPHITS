from collections import namedtuple, Iterable

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



class PhitsObject:
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

    There's also a series of parameters related to grouping, that allow
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
               "group_by", "separator", "prelude", "max_groups"}
    names = {"parameters", "source", "material", "surface", "cell", "transform", "temperature","mat_time_change","magnetic_field",
             "neutron_magnetic_field", "mapped_magnetic_field", "uniform_electromagnetic_field", "mapped_electromagnetic_field",
             "delta_ray", "track_structure", "super_mirror", "elastic_option", "importance", "weight_window", "ww_bias",
             "forced_collisions", "repeated_collisions", "volume", "multiplier", "mat_name_color", "reg_name", "counter", "timer",
             "t-track", "t-cross", "t-point", "t-adjoint", "t-deposit", "t-deposit2", "t-heat", "t-yield", "t-product", "t-dpa",
             "t-let", "t-sed", "t-time", "t-interact", "t-dchain", "t-wwg", "t-wwbg", "t-volume", "t-gshow", "t-rshow","t-3dshow"}

    def __init__(self, *args,  **kwargs):
        assert self.name in self.names, f"Unrecognized PHITS type {self.name} in PhitsObject initialization."

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

    def prelude(self):
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
    def __init__(self, **kwargs):
        self.name = "parameters"
        for k,v in kwargs.items():
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

class Mesh():
    """Represents all list-typed data in PHITS."""
    def __init__(self, axis, bins=None): # All bin generation is easily done in Python via list comprehensions
        assert axis in ["energy", "time", "x", "y", "z", "radius", "angle", "let"], f"Unrecognized axis {axis} in mesh definition."
        self.axis = axis
        self.bins = tuple(bins)
        print(self.bins)
        if axis != "angle":
            self.type = 2
        else:
            pass  # TODO: figure out angle mesh


    def __eq__(self, other):
        if type(self) != type(other):
            return False

        else:
            return {k: v for k, v in self.__dict__.items() if v is not other} \
                == {k: v for k, v in other.__dict__.items() if v is not self}

    def __hash__(self):
        return hash(tuple(v for k, v in sorted(self.__dict__.items()) \
                          if (self not in v.__dict__.values() if hasattr(v, "__dict__") else True)))
