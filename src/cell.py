from copy import deepcopy
from base import *
from transform import *
from misc import *
from material import *
from surface import TetrahedronBox, surface_spec



subobject_syntax = {"magnetic_field": (None, OneOf(IsA(MagneticField, index=True), IsA(NeutronMagneticField, index=True), #), None),
                                                   IsA(MappedMagneticField, index=True)), None),
                    "electromagnetic_field": (None, OneOf(IsA(UniformElectromagneticField, index=True),#), None),
                                                          IsA(MappedElectromagneticField, index=True)), None),
                    "delta_ray": (None, IsA(DeltaRay, index=True), None),
                    "track_structure": (None, IsA(TrackStructure, index=True), None),
                    # "super_mirror": (None, IsA(SuperMirror, index=True), None),
                    "elastic_option": (None, IsA(ElasticOption, index=True), None),
                    "importance": (None, IsA(Importance, index=True), None),
                    "weight_window": (None, IsA(WeightWindow, index=True), None),
                    "ww_bias": (None, IsA(WWBias, index=True), None),
                    "forced_collisions": (None, IsA(ForcedCollisions, index=True), None),
                    "repeated_collisions": (None, IsA(RepeatedCollisions, index=True), None),
                    "reg_name": (None, IsA(RegionName, index=True), None),
                    "counter": (None, IsA(Counter, index=True), None),
                    "timer": (None, IsA(Timer, index=True), None)}

common_syntax = subobject_syntax | {"volume": ("VOL", PosReal(), None),
                                    "temperature": ("TMP", PosReal(), None),
                                    "transform": ("TRCL", IsA(Transform, index=True), None)}

class Tetrahedral(PhitsObject):
    """A box filled with tetrahedrons."""
    name = "cell"
    syntax = common_syntax | {"regions": (None, IsA(TetrahedronBox, index=True), 0),
                              "material": (None, IsA(Material, index=True), 1),
                              "density": (None, PosReal(), 2),
                              "tet_format": (None, FinBij({"tetgen": "tetgen", "NASTRAN": "NASTRAN"}), 1),
                              "tet_file": (None, Path(), 2),
                              "scale_factor": ("TSFAC", PosReal(), None)}

    shape = lambda self: (("self", "material", "density", "regions", "\\"),
                          "volume\\", "temperature\\", "transform\\", "LAT=3\\",
                          f"tfile={self.tet_file}" if self.tet_format == "tetgen" else f"nfile={self.tet_file}", "scale_factor")

    subobjects = set(subobject_syntax.keys())

    # def restrictions(self):
    #     if len(self.regions) != 1:
    #         raise ValueError(f"Tetrahedral cells may have only one TetrahedronBox region; got {self.regions}")
        # if self.forced_collisions is not None and self.repeated_collisions is not None:
        #     raise ValueError(f"Cannot set both forced_collisions and repeated_collisions on a Tetrahedral cell.")

    def __or__(self, other): # Union of cells; adopts leftmost's properties
        r = deepcopy(self)
        setattr(r, "regions", self.regions + ("|",) + other.regions)
        return r

    def __invert__(self): # Set complement of cell; new cell has old properties
        r = deepcopy(self)
        r.regions = ("~", self.regions)
        return r

    def __and__(self, other): # Intersection of cells; drops properties
        r = deepcopy(self)
        r.regions = self.regions + other.regions
        return r

    def __rshift__(self, other): # returns other's regions with self's properties
        r = deepcopy(self)
        r.regions = other.regions
        return r

    def __lshift__(self, other): # returns self's region with other's properties
        r = deepcopy(other)
        r.regions = self.regions
        return r



class Void(PhitsObject):
    """A region with no material, just vacuum."""
    name = "cell"
    syntax = common_syntax | {"regions": (None, RegionTuple(surface_spec), 0)}
    shape = lambda self: (("self", "0", "regions", "\\"), "volume\\", "temperature\\", "transform\\", "")
    subobjects = set(subobject_syntax.keys())

    def __or__(self, other): # Union of cells; adopts leftmost's properties
        r = deepcopy(self)
        setattr(r, "regions", (self.regions,) + ("|",) + (other.regions,))
        return r

    def __invert__(self): # Set complement of cell; new cell has old properties
        r = deepcopy(self)
        r.regions = ("~", (self.regions,))
        return r

    def __and__(self, other): # Intersection of cells; drops properties
        r = deepcopy(self)
        r.regions = (self.regions,) + (other.regions,)
        return r

    def __rshift__(self, other): # returns other's regions with self's properties
        r = deepcopy(self)
        r.regions = other.regions
        return r

    def __lshift__(self, other): # returns self's region with other's properties
        r = deepcopy(other)
        r.regions = self.regions
        return r

    # def restrictions(self):
    #     if len(self.regions) != 1:
    #         raise ValueError(f"Tetrahedral cells may have only one TetrahedronBox region; got {self.regions}")
        # if self.forced_collisions is not None and self.repeated_collisions is not None:
        #     raise ValueError(f"Cannot set both forced_collisions and repeated_collisions on a Tetrahedral cell.")


class OuterVoid(PhitsObject):
    """Void, but different for some reason. Probably shouldn't be used directly;
    `run_phits` creates the required OuterVoid for you automatically."""
    name = "cell"
    syntax = common_syntax | {"regions": (None, RegionTuple(surface_spec), 0)}
    shape = lambda self: (("self", "-1", "regions", "\\"), "volume\\", "temperature\\", "transform\\", "")
    subobjects = set(subobject_syntax.keys())

    def __or__(self, other): # Union of cells; adopts leftmost's properties
        r = deepcopy(self)
        setattr(r, "regions", (self.regions,) + ("|",) + (other.regions,))
        return r

    def __invert__(self): # Set complement of cell; new cell has old properties
        r = deepcopy(self)
        r.regions = ("~", (self.regions,))
        return r

    def __and__(self, other): # Intersection of cells; drops properties
        r = deepcopy(self)
        r.regions = (self.regions,) + (other.regions,)
        return r
    def __rshift__(self, other): # returns other's regions with self's properties
        r = deepcopy(self)
        r.regions = other.regions
        return r

    def __lshift__(self, other): # returns self's region with other's properties
        r = deepcopy(other)
        r.regions = self.regions
        return r

# TODO: operations
class Cell(PhitsObject):
    """The prototypical `Cell`, consisting of the intersection of several regions defined by surfaces with a material of some density."""
    name = "cell"
    syntax = common_syntax | {"regions": (None, RegionTuple(surface_spec), 0),
                              "material": (None, IsA(Material, index=True), 1),
                              "density": (None, PosReal(), 2)}
    shape = lambda self: (("self", "material", "density", "regions", "\\"),
                          "volume\\", "temperature\\", "transform\\", "")

    subobjects = set(subobject_syntax.keys())

    def __or__(self, other): # Union of cells; adopts leftmost's properties
        r = deepcopy(self)
        setattr(r, "regions", (self.regions,) + ("|",) + (other.regions,))
        return r

    def __invert__(self): # Set complement of cell; new cell has old properties
        r = deepcopy(self)
        r.regions = ("~", (self.regions,))
        return r

    def __and__(self, other): # Intersection of cells; drops properties
        r = deepcopy(self)
        r.regions = (self.regions,) + (other.regions,)
        return r

    def __rshift__(self, other): # returns other's regions with self's properties
        r = deepcopy(self)
        r.regions = other.regions
        return r

    def __lshift__(self, other): # returns self's region with other's properties
        r = deepcopy(other)
        r.regions = self.regions
        return r

# idea: generate a UUID for the universe/fill, and then map UUIDs -> index at runtime
# other idea: make a Universe class, define an __init__, and make a call to super() for the normal __init__,
# but use the rest of __init__ to set the right attributes on the underlying cells, and marshall definitions
# def fill_universe(mask: Cell, contents: list[Cell]):
#     pass
