from valspec import *
import re

__pdoc__ = dict()
__pdoc__["builds"] = False
__pdoc__["slices"] = False

test = 1
"""
do   | tables
- | ---
even | work
"""

def _continue_lines(inp: str) -> str:
    """If a line is too long for PHITS to handle, use PHITS's line continuation syntax to fix it."""
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
readable_remapping = True
"""If `True`, use a remapped initialization syntax that's more informative.
   E.g., `Parameters(control="output_echo_only")` instead of `Parameters(icntl=3)`.

   If `False`, use a syntax as close as possible to that specified in the PHITS manual (i.e., the latter form in the above example).
   Note that some identifiers in PHITS do not conform to Python's identifier syntax (e.g. 2d-type, emin(14));
   these identifiers are sanitized as follows:
       - dashes -> underscores
       - beginning with a number -> that clause moved to the end
       - parentheses -> omitted"""


class PhitsObject:
    """The base factory class distinguishing objects that are intended to end up in some section of a `.inp` file,
    and defining equality and hashability of such objects in sensible ways.

    PhitsObject values always correspond to some section of an input file.
    """

    name = None
    """A string corresponding to the PHITS .inp section an object appears in. See `PhitsObject.names`"""

    syntax = dict()
    """A dictionary with entries of the form `"python_identifier": ("PHITS_identifier", acceptable_ValSpec, arg_position, Optional(none_val))`.
    The key is the attribute on the Python PhitsObject instance, and the keyword argument necessary to set it.
    The second is what type the object must be, as a ValSpec.
    The third is what index in `*args` the argument must have; if it is set to `None`, the argument is optional.
    The last is a value to put in the .inp if the attribute is `None` at compile-time; ordinarily, it's just nothing.
    The first two arguments can be tuples, in which case the passed value must be a tuple of the specified type;
    the single Python assignment corresponds to the entrywise assignments of this tuple to the PHITS identifiers in the .inp."""

    shape = tuple()
    """A tuple that details how the object is to be represented in the .inp file.
    The syntax is inspired by Emacs Lisp skeletons, with some optimizations to tailor it for this use-case.
    In general, the strings in the tuple are inserted verbatim, with a newline in between them.
    If an entry is the name of an attribute of the instance, say `"attr"` whose value is not `None`,
    then the string `PHITS_identifier = acceptable_ValSpec.phits(self.attr)`, based on a lookup of `attr` in `PhitsObject.syntax`, is inserted.
    If an entry is `self`, then `self.index` is inserted.
    To avoid the above two behaviors, and insert a string verbatim, prepend a quote `'`.
    If an entry is another tuple, that tuple gets evaluated as above, but the `identifier = ` is not inserted,
    and a mere space separates the entries.
    Appending a `\\` to any string disables the insertion of the spacing that would otherwise follow;
    similarly, having `\\` as the last string of a tuple entry disables the newline that would otherwise follow.

    If this attribute is callable, it gets called on the instance of the PhitsObject, and the result is processed according to the rules above.
    """

    index = None
    """The compile-time assigned number of an object. Should not be set directly."""

    no_hash = {"index", "value_map", "ident_map", "nones", "shape", "subobjects", "required", "positional", "optional",
               "group_by", "separator", "prelude", "max_groups", "group_size", "parser", "validator"}
    """Attributes that don't affect the identity of a PhitsObject."""

    names = {"parameters", "source", "material", "surface", "cell", "transform", "temperature","mat_time_change","magnetic_field",
             "neutron_magnetic_field", "mapped_magnetic_field", "uniform_electromagnetic_field", "mapped_electromagnetic_field",
             "delta_ray", "track_structure", "super_mirror", "elastic_option", "importance", "weight_window", "ww_bias",
             "forced_collisions", "repeated_collisions", "volume", "multiplier", "mat_name_color", "reg_name", "counter", "timer",
             "t-track", "t-cross", "t-point", "t-adjoint", "t-deposit", "t-deposit2", "t-heat", "t-yield", "t-product", "t-dpa",
             "t-let", "t-sed", "t-time", "t-interact", "t-dchain", "t-wwg", "t-wwbg", "t-volume", "t-gshow", "t-rshow","t-3dshow"}
    """The permissible `PhitsObject.name`s."""

    group_by = None
    """A key function by which to group PhitsObjects together.
    These groups are necessary, for instance, in the `Importance` section, where one can set different importances for different particles
    in the same cell."""

    separator = None
    """A function returning a string to be placed between the definitions of the groups; see `PhitsObject.group_by`."""

    max_groups = None
    """The maximum number of groups of objects of a given type; see `PhitsObject.group_by`."""

    prelude = tuple()
    """A skeleton just like `PhitsObject.shape`, but inserted before all definitions of the `PhitsObject` subclass in question."""
    def __init__(self, *args,  **kwargs):
        """Arguments are interpreted according to `PhitsObject.syntax`, and then any leftovers in `kwargs` are used to create a `Parameters`
        object, which is then assigned to a `parameters` attribute."""
        assert self.name in self.names, f"Unrecognized PHITS type {self.name} in PhitsObject initialization."
        required = list(map(lambda tup: tup[0],
                            sorted([(k, v) for k, v in self.syntax.items() if v[2] is not None], key=lambda tup: tup[1][2])))
        assert len(args) == len(required), f"Wrong number of positional arguments specified in the definition of {self.name} object."
        for idx, arg in enumerate(args):
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



    def _add_definition(self, how: tuple, to: str, assignments: bool = True) -> str:
        """Recursively performs skeleton insertions according to `how` at the end of `to`."""
        if callable(how):
            how = how(self)

        for attr in how:
            if isinstance(attr, str) and attr[0] == "'":
                to += attr[1:]
                return to

            endstr = "\n" if assignments else " "
            spacing = " "

            if isinstance(attr, str) and attr[-1] == "\\":
                endstr = " "
                spacing = ""
                attr = attr[:-1]

            if isinstance(attr, tuple):
                to += self.add_definition(attr, to, assignments=False)
                if attr[-1] == "\\":
                    to += " "
                else:
                    to += "\n"

            elif attr in self.syntax:
                val = getattr(self, attr)
                phits_iden = self.syntax[attr][0]
                valspec = self.syntax[attr][1]
                noneval = ""
                if len(self.synatx[attr]) == 4:
                    noneval = self.syntax[attr][3]

                if val is not None:
                    if isinstance(phits_iden, tuple):
                        for i, (phits, spec) in enumerate(zip(phits_iden, valspec)):
                            assign = f"{phits}{spacing}={spacing}" if assignments else ""
                            to += f"{assign}{spec.phits(val[i])}{endstr}"
                    else:
                        assign = f"{phits_iden}{spacing}={spacing}" if assignments else ""
                        to += f"{assign}{valspec.phits(val)}{endstr}"

            else:
                if attr == "self":
                    to += self.index
                    to += endstr
                else:
                    to += attr
                    to += endstr
        return to

    def prelude_str(self) -> str:
        """Return a string to appear before the collection of all definitions of subclass instances in an `.inp` file."""
        inp = self.add_definition(self.prelude, "")

        return _continue_lines(inp)

    def definition(self) -> str:
        """Return the string representing the particular PhitsObject in an `.inp` file."""
        inp = self.add_definition(self.shape, "")

        return _continue_lines(inp)

    def section_title(self) -> str:
        """Return the section title under which a PhitsObject belongs."""
        sec_name = self.name.replace("_", " ").title()
        return f"[{sec_name}]\n"

    @classmethod
    def syntax_desc(self) -> str:
        """Return a readable summary of the initialization syntax of the PhitsObject in question.
        Used to generate documentation, but is useful in interactive sessions to """
        required = sorted([(k, v) for k, v in self.syntax.items() if v[2] is not None], key=lambda tup: tup[1][2])
        opt = [(k, v) for k, v in self.syntax.items() if v[2] is None]
        r = ""
        if required:
            r = "Required arguments:\nPHITS name|Accepted value\n----|----\n"
            for py_attr, (phits_attr, valspec, position) in required:
                r += f"{phits_attr}|{valspec.description()}\n"

        if opt:
            r += "Optional arguments:\nPython name|PHITS name|Accepted value\n----|----|----\n"
            for py_attr, (phits_attr, valspec, position) in opt:
                if isinstance(valspec, tuple):
                    j = ", "
                    r += f"`{py_attr}`|`{phits_attr}`|A tuple ({j.join(map(lambda x: x.description().capitalize(), valspec))}).\n"
                else:
                    r += f"`{py_attr}`|`{phits_attr}`|{valspec.description().capitalize()}.\n"
        return r


    @classmethod
    def syntax_for(self, attr: str) -> str:
        """Return a readable summary of a specific initialization parameter of the PhitsObject in question."""
        r = "PHITS name\tAccepted value\tPosition\n"
        r += f"{self.syntax[attr][0]}\t{self.syntax[attr][1]}\t{self.syntax[attr][2]}"
        return r

    def _sanitize(self, iden: str) -> str:
        """Make a PHITS identifier Python-acceptable.
        See `readable_syntax`.
        """
        r = iden.replace("-", "_")
        r = r.replace("(", "")
        r = r.replace(")", "")
        if r[0].isdigit():
            if m := re.match(".*?_", r):
                r = r[m.end():] + "_" + r[:m.end() - 1]

        return r


    def __eq__(self, other: "PhitsObject") -> bool:
        """PhitsObjects should be equal if their definitions would be the same."""
        if type(self) != type(other):
            return False
        elif hasattr(self, "__dict__") and hasattr(other, "__dict__"):
            d1 = {k: v for k, v in self.__dict__.items() if k not in self.no_hash}
            d2 = {k: v for k, v in other.__dict__.items() if k not in self.no_hash}
            return d1 == d2

    def __hash__(self) -> int:
        """PhitsObjects have the hash of their identity-defining attributes."""
        return hash(tuple(v for k, v in sorted(self.__dict__.items()) \
                          if (self not in v.__dict__.values() if hasattr(v, "__dict__") else True) and k not in self.no_hash))
