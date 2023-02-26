# I want a way to parse and unparse.

# E.g., I have a string representation of a Python object, and a parser/transformer combo that translates that into said object.
# Suppose I want to go the other way---given a Python object, produce a string according to the grammar that represents it.

# The Python object should be given information about the abstract definitional properties of the object (i.e. the attributes),
# and generate the ways in which it's constructed from Python and from a string.
# There should be support for mapping the attribute names in different ways when creating and string-representing,
# and so also for the attribute values.

# In each direction, you need a consumption syntax and a production syntax, specifying:
#    - The way

# current idea: define an __init__ which produces and sets a parse-tree representation of the object;
# one can use Lark's Reconstructor to produce the text.

class Bijection(Transformer, Parser):

class PhitsObject(type):

class Void(PhitsObject):

    attributes = {"name": "cell",}
    name = "cell"
    required = ["regions"]
    positional = ["regions"]
    optional = ["transform", "temperature", "magnetic_field", "neutron_magnetic_field",
                "mapped_magnetic_field", "uniform_electromagnetic_field", "mapped_electromagnetic_field",
                "delta_ray", "track_structure", "super_mirror", "elastic_option", "importance",
                "weight_window", "ww_bias", "forced_collisions", "repeated_collisions", "volume",
                "reg_name", "counter", "timer", "tally", "containing_universe", "lattice", "universe_contents",
                "tet_format", "tet_file", "tet_scale"]
    grammar = r"""

    """
    shape = (("self", "-1\\"), lambda self: tup_to_def(self.regions),
             "volume\\", "temperature\\", "transform\\", "containing_universe\\", "lattice\\",
             lambda self: ("LAT=3 " + (f"tfile={self.tet_file}" if self.tet_format == "tetgen" else f"nfile={self.tet_file}")) if self.lattice is not None else "",
             "tet_scale\\", lambda self: f"FILL={self.index}" if self.universe_contents else "")
    ident_map = {"volume": "VOL", "temperature": "TMP", "transform": "TRCL", "containing_universe": "U", "lattice": "LAT",
                 "tet_scale": "TSFAC"},
    value_map = {"rectangular": 1, "hexagonal": 2}
    subobjects = ["transform", "temperature", "magnetic_field", "neutron_magnetic_field",
                  "mapped_magnetic_field", "uniform_electromagnetic_field", "mapped_electromagnetic_field",
                  "delta_ray", "track_structure", "super_mirror", "elastic_option", "importance",
                  "weight_window", "ww_bias", "forced_collisions", "repeated_collisions", "volume",
                  "reg_name", "counter", "timer", "tally"]
    parser = r"""
    start =
"""

