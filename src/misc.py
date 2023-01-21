from base import *

class MagneticField(PhitsObject): # Right now, the only way to set this is to do setattr(<cell>, "magnetic_field", <MagneticField>) since one cannot pass the cell to create the magnetic field while the cell is being initialized.
    name = "magnetic_field"
    required = ["typ", "gap", "strength"]
    positional = ["typ", "gap", "strength"]
    optional = ["transform", "time", "cell"]
    shape = (("cell", "typ", "gap", "strength", "transform", "time"))
    prelude = (("reg", "\\typ", "\\gap", "mgf", "trcl", "\\time"))
    nones = {"transform": 0, "time": "non"}

class NeutronMagneticField(PhitsObject):
    name = "neutron_magnetic_field"
    required = ["typ", "gap", "strength"]
    positional = ["typ", "gap", "strength"]
    optional = ["transform", "polar", "time", "cell"]
    shape = (("cell", "typ", "gap", "strength", "transform", "polar", "time"))
    prelude = (("reg", "\\typ", "\\gap", "mgf", "trcl", "\\polar", "time"))
    nones = {"transform": 0, "polar": "non", "time": "non"}

class MappedMagneticField(PhitsObject):
    name = "neutron_magnetic_field"
    required = ["typ", "gap", "strength", "m_file"]
    positional = ["typ", "gap", "strength", "m_file"]
    optional = ["transform", "cell"]
    shape = (("cell", "typ", "gap", "strength", "transform", "m_file"))
    prelude = (("reg", "\\typ", "\\gap", "mgf", "trcl", "file"))
    nones = {"transform": 0}


class UniformElectromagneticField(PhitsObject):
    name = "uniform_electromagnetic_field"
    required = ["e_strength", "m_strength"]
    positional = ["e_strength", "m_strength"]
    optional = ["e_transform", "m_transform"]
    shape = (("cell", "e_strength", "m_strength", "e_transform", "m_transform"))
    prelude = (("reg", "elf", "mgf", "trcle", "trclm"))
    nones = {"e_transform": 0, "m_transform": 0}


class MappedElectromagneticField(PhitsObject):
    name = "mapped_electromagnetic_field"
    required = ["typ_e", "typ_m", "gap", "e_strength", "m_strength", "e_file", "m_file"]
    positional = ["typ_e", "typ_m", "gap", "e_strength", "m_strength", "e_file", "m_file"]
    optional = ["e_transform", "m_transform"]
    shape = (("cell", "typ_e", "typ_m", "gap", "e_strength", "m_strength", "e_transform", "m_transform", "e_file", "m_file"))
    prelude = (("reg", "type", "typm", "gap", "elf", "mgf", "trcle", "trclm", "filee", "filem"))
    nones = {"e_transform": 0, "m_transform": 0}


class DeltaRay(PhitsObject):
    name = "delta_ray"
    required = ["threshold"]
    positional = ["threshold"]
    optional = ["cell"]
    shape = (("cell", "threshold"))
    prelude = (("reg", "del"))


class TrackStructure(PhitsObject):
    name = "track_structure"
    required = ["mID"]
    positional = ["mID"]
    optional = ["cell"]
    shape = (("cell", "mID"))
    prelude = ("reg", "\\mID")


class SuperMirror(PhitsObject):
    name = "super_mirror"
    required = ["r_in", "r_out", "mirror_material", "reflectivity", "cutoff", "falloff_rate", "cutoff_width"]
    positional = ["r_in", "r_out", "mirror_material", "reflectivity", "cutoff","falloff_rate", "cutoff_width"]
    optional = []
    shape = (("r_in", "r_out", "mirror_material", "reflectivity", "cutoff", "falloff_rate", "cutoff_width"))
    prelude = (("r-in", "r-out", "mm", "r0", "qc", "am", "wm"))


class ElasticOption(PhitsObject):
    name = "elastic_option"
    required = ["c1", "c2", "c3", "c4"]
    positional = ["c1", "c2", "c3", "c4"]
    optional = ["cell"]
    shape = (("cell", "c1", "c2", "c3", "c4"))
    prelude = (("reg", "\\c1", "\\c2", "\\c3", "\\c4"))
    nones = {"c1": "non", "c2": "non", "c3": "non", "c4": "non"}


class FragData(PhitsObject):
    name = "frag_data"
    required = ["option", "projectile", "target", "file"]
    positional = ["option", "projectile", "target", "file"]
    optional = []
    shape = (("option", "projectile", "target", "file"))
    prelude = (("opt", "proj", "targ", "\\file"))


class Importance(PhitsObject):
    name = "importance"
    required = ["particles", "importance"]
    positional = ["particles", "importance"]
    optional = ["cell"]
    shape = (("cell", "importance"))
    prelude = ("particles", ("reg", "imp"))
    group_by = lambda self: self.particles
    separator = lambda self: self.section_title()
    max_groups = 6
    ident_map = {"particles": "part"}


class WeightWindow(PhitsObject):
    name = "weight_window"
    required = ["mesh", "particles", "windows"]
    positional = ["mesh", "particles", "windows"]
    optional = ["grid", "cell"]
    shape = (("cell", "windows"))
    prelude = ("mesh", "particles", "grid", ("reg", lambda self: " ".join(f"ww{i}" for i in range(1, len(self.windows) + 1))))
    group_by = lambda self: (self.mesh, self.particles, self.grid)
    separator = lambda self: self.section_title()
    max_groups = 6
    ident_map = {"particles": "part"}


class WWBias(PhitsObject):
    name = "ww_bias"
    required = ["particles", "mesh", "biases"]
    positional = ["particles", "mesh", "biases"]
    optional = ["cell"]
    shape = (("cell", "biases"))
    prelude = ("particles", "mesh", ("reg", lambda self: " ".join(f"wwb{i}" for i in range(1, len(self.biases) + 1))))
    group_by = lambda self: (self.particles, self.mesh)
    separator = lambda self: self.section_title()
    max_groups = 6
    ident_map = {"particles": "part"}


class ForcedCollisions(PhitsObject):
    name = "forced_collisions"
    required = ["particles", "probability"]
    positional = ["particles", "probability"]
    optional = ["cell"]
    shape = (("cell", "probability"))
    prelude = ("particles", ("reg", "fcl"))
    group_by = lambda self: self.particles
    separator = lambda self: self.section_title()
    max_groups = 6
    ident_map = {"particles": "part"}



class RepeatedCollisions(PhitsObject):
    name = "repeated_collisions"
    required = ["particles", "mother", "ebounds", "collision_reps", "evaporation_reps"]
    positional = ["particles", "mother", "ebounds", "collision_reps", "evaporation_reps"]
    optional = ["cell"]
    shape = (("cell", "collision_reps", "evaporation_reps"))
    prelude = ("particles", "mother", "ebounds", ("reg", "n-coll", "n-evaps"))
    group_by = lambda self: (self.particles, self.mother)
    separator = lambda self: self.section_title()
    max_groups = 6
    ident_map = {"particles": "part", "ebounds": ("emin", "emax")}


class Volume(PhitsObject):
    name = "volume"
    required = ["volume"]
    positional = ["volume"]
    optional = ["cell"]
    shape = (("cell", "volume"))
    prelude = (("reg", "vol"))


# TODO: finish this
class Multiplier(PhitsObject):
    name = "multiplier"
    required = ["number", "interpolation", "particles"]
    positional = ["number", "interpolation", "particles"]
    optional = []



class RegionName(PhitsObject):
    name = "region_name"
    required = ["reg_name", "size"]
    positional = ["reg_name", "size"]
    optional = ["cell"]
    shape = (("cell", "reg_name", "size"))



class Counter(PhitsObject):
    name = "counter"
    required = ["particles", "in", "out", "collisions", "reflections"]
    positional = ["particles", "in", "out", "collisions", "reflections"]
    optional = ["cell", "particles_except"]
    shape = (("cell", "in", "out", "collisions", "reflections"))
    prelude = ("particles", "particles_except", ("reg", "in", "out", "coll", "ref"))
    group_by = lambda self: self.particles
    separator = lambda self: f"counter = {self.index}"
    ident_map = {"particles": "part", "particles_except": "*part"}
    max_groups = 3

class Timer(PhitsObject):
    name = "timer"
    required = ["enter", "out", "coll", "ref"]
    positional = ["enter", "out", "coll", "ref"]
    optional = ["cell"]
    shape = (("cell", "enter", "out", "coll", "ref"))
    prelude = (("reg", "in", "out", "coll", "ref"))
