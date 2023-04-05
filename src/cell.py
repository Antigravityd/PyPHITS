import copy
from base import *


def tup_to_def(tup):
    r = ""
    for el in tup:
        if isinstance(el, PhitsObject):
            sense = "-" if el.inside is not None else ""
            r += f"{sense}{el.index} "
        elif isinstance(el, tuple):
            r += "(" + tup_to_def(el) + ") "
        elif el == "~":
            r += "#"
        elif el == "|":
            r += ": "
        else:
            breakpoint()
            raise ValueError(f"Unrecognized token {el} in cell region definition.")

    return r


# TODO: define universe object
# "containing_universe": ("U", PosInt(), None), # should be dynamically assigned
subobject_syntax = {"universe": ("FILL", List(IsA(Cell, index=True)), None),
                    "transform": ("TRCL", IsA(Transform, index=True), None),
                    "temperature": (None, IsA(Temperature), None),
                    "magnetic_field": (None, IsA(MagneticField), None),
                    "neutron_magnetic_field": (None, IsA(NeutronMagneticField), None),
                    "mapped_magnetic_field": (None, IsA(MappedMagneticField), None),
                    "uniform_electromagnetic_field": (None, IsA(UniformElectromagneticField), None),
                    "mapped_electromagnetic_field": (None, IsA(MappedElectromagneticField), None),
                    "delta_ray": (None, IsA(DeltaRay), None),
                    "track_structure": (None, IsA(TrackStructure), None),
                    "super_mirror": (None, IsA(SuperMirror), None),
                    "elastic_option": (None, IsA(ElasticOption), None),
                    "importance": (None, IsA(Importance), None),
                    "weight_window": (None, IsA(WeightWindow), None),
                    "ww_bias": (None, IsA(WWBias), None),
                    "forced_collisions": (None, IsA(ForcedCollisions), None),
                    "repeated_collisions": (None, IsA(RepeatedCollisions), None),
                    "volume": (None, IsA(Volume), None),
                    "reg_name": (None, IsA(RegName), None),
                    "counter": (None, IsA(Counter), None),
                    "timer": (None, IsA(Timer), None),
                    "tally": (None, IsA(Tally), None)}

common_syntax = subobject_syntax | {"containing_universe": ("U", PosInt(), None),
                                    "lattice": ("LAT", FinBij({"quadrilateral": 1, "hexagonal": 2, "tetrahedral": 3}), None),
                                    "tet_format": (None, )}

class Tetrahedral(PhitsObject):
    name = "cell"
    syntax = common_syntax | {"within": (None, IsA(TetrahedronBox, index=True), 0),
                              "tet_format": (None, FinBij({"tetgen": "tetgen", "NASTRAN": "NASTRAN"}), 1),
                              "tet_file": (None, Path(), 2),
                              "scale_factor": ("TSFAC", PosReal(), None)}

    shape = lambda self: (("self", "material", "density\\"), tup_to_def((self.within)),
                          "volume\\", "temperature\\", "transform\\", f"U={self.universe}" if self.universe is not None else "",
                          "LAT=3", f"tfile={self.tet_file}" if self.tet_format == "tetgen" else f"nfile={self.tet_file}",
                          "scale_factor\\")

class Lattice(PhitsObject):
    pass


class Void(PhitsObject):
    name = "cell"
    syntax = {"regions": (None, List(Tuple(IsA(Surface, index=True), Orientation())), 0),
              "containing_universe": ("U", IsA(Cell, index=True), None),
              "lattice": ("LAT", ), "universe_contents": tuple()}
    required = ["regions"]
    positional = ["regions"]
    optional = ["transform", "temperature", "magnetic_field", "neutron_magnetic_field",
                "mapped_magnetic_field", "uniform_electromagnetic_field", "mapped_electromagnetic_field",
                "delta_ray", "track_structure", "super_mirror", "elastic_option", "importance",
                "weight_window", "ww_bias", "forced_collisions", "repeated_collisions", "volume",
                "reg_name", "counter", "timer", "tally", "containing_universe", "lattice", "universe_contents",
                "tet_format", "tet_file", "tet_scale"]
    shape = (("self", "-1\\"), lambda self: tup_to_def(self.regions),
             "volume\\", "temperature\\", "transform\\", "containing_universe\\", "lattice\\",
             lambda self: ("LAT=3 " + (f"tfile={self.tet_file}" if self.tet_format == "tetgen" else f"nfile={self.tet_file}")) if self.lattice is not None else "",
             "tet_scale\\", lambda self: f"FILL={self.index}" if self.universe_contents else "")
    ident_map = {"volume": "VOL", "temperature": "TMP", "transform": "TRCL", "containing_universe": "U", "lattice": "LAT",
                 "tet_scale": "TSFAC"},
    value_map = {"rectangular": 1, "hexagonal": 2}
    subobjects = set(subobject_syntax.keys())




# the regions should be a tuple of surfaces with the logical operators interspersed as strings. Sub-tuples indicate parentheses,
# and the surfaces are interpreted in the "positive" sense just like in PHITS
class Cell(PhitsObject):
    name = "cell"
    required = ["regions", "material", "density"]
    positional = ["regions", "material", "density"]
    optional = ["transform", "temperature", "magnetic_field", "neutron_magnetic_field",
                "mapped_magnetic_field", "uniform_electromagnetic_field", "mapped_electromagnetic_field",
                "delta_ray", "track_structure", "super_mirror", "elastic_option", "importance",
                "weight_window", "ww_bias", "forced_collisions", "repeated_collisions", "volume",
                "reg_name", "counter", "timer", "tally", "containing_universe", "lattice", "universe_contents",
                "tet_format", "tet_file", "tet_scale"]
    shape = (("self", "material", "density\\"), lambda self: tup_to_def(self.regions),
             "volume\\", "temperature\\", "transform\\", "containing_universe\\", "lattice\\",
             lambda self: "LAT=3 " + (f"tfile={self.tet_file}" if self.tet_format == "tetgen" else f"nfile={self.tet_file}") \
             if self.lattice is not None else "",
             "tet_scale\\", lambda self: f"FILL={self.index}" if self.universe_contents else "")
    ident_map = {"volume": "VOL", "temperature": "TMP", "transform": "TRCL", "containing_universe": "U", "lattice": "LAT",
                 "tet_scale": "TSFAC"},
    value_map = {"rectangular": 1, "hexagonal": 2}
    subobjects = ["transform", "temperature", "magnetic_field", "neutron_magnetic_field",
                  "mapped_magnetic_field", "uniform_electromagnetic_field", "mapped_electromagnetic_field",
                  "delta_ray", "track_structure", "super_mirror", "elastic_option", "importance",
                  "weight_window", "ww_bias", "forced_collisions", "repeated_collisions", "volume",
                  "reg_name", "counter", "timer", "tally"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.universe_contents:
            for cell in self.fill:
                cell.containing_universe = self


    def __or__(self, other): # Union of cells; adopts leftmost's properties
        r = copy.deepcopy(self)
        setattr(r, "regions", self.regions + ("|",) + other.regions)
        return r

    def __inv__(self): # Set complement of cell; new cell has old properties
        r = copy.deepcopy(self)
        r.regions = ("~", self.regions)
        return r

    def __and__(self, other): # Intersection of cells; drops properties
        r = copy.deepcopy(self)
        r.regions = self.regions + other.regions
        return r

    def __rshift__(self, other): # returns other's regions with self's properties
        r = copy.deepcopy(self)
        r.regions = other.regions
        return r

    def __lshift__(self, other): # returns self's region with other's properties
        r = copy.deepcopy(other)
        r.regions = self.regions
        return r
