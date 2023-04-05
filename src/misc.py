from base import *

# no temperature; do that at the cell level for now

class MagneticField(PhitsObject): # Right now, the only way to set this is to do setattr(<cell>, "magnetic_field", <MagneticField>) since one cannot pass the cell to create the magnetic field while the cell is being initialized.
    name = "magnetic_field"
    syntax = {"typ": (None, FinBij({"dipole": 2, "quadrupole": 4}), 0),
              "strength": (None, Real(), 1),
              "gap": (None, PosReal(), None, 0.0),
              "transform": (None, IsA(Transform, index=True), None, 0),
              "time": (None, PosReal(), None, "non"),
              "cell": (None, IsA(Cell, index=True), None)}
    prelude = (("reg", "'typ", "'gap", "mgf", "trcl", "'time"))
    shape = (("cell", "typ", "gap", "strength", "transform", "time"))


class NeutronMagneticField(PhitsObject):
    name = "neutron_magnetic_field"
    syntax = {"typ": (None, FinBij({"identified": 60, "nograv": 61, "dipole": 102,
                                    "quadrupole": 104, "sextupole": 106}), 0),
              "strength": (None, Real(), 1),
              "gap": (None, PosReal(), None, 0.0),
              "polarization": (None, Real(), "non")
              "transform": (None, IsA(Transform, index=True), None, 0),
              "time": (None, PosReal(), None, "non"),
              "cell": (None, IsA(Cell, index=True), None)}
    prelude = (("reg", "'typ", "'gap", "mgf", "trcl", "'polar", "time"))
    shape = (("cell", "typ", "gap", "strength", "transform", "polar", "time"))

class MappedMagneticField(PhitsObject):
    name = "neutron_magnetic_field"
    syntax = {"typ": (None, FinBij({"xyz_list_charged": -1, "rz_list_charged": -2, "xyz_map_charged": -3, "rz_map_charged": -4,
                                    "xyz_list_neutron": -101, "rz_list_neutron": -102, "xyz_map_neutron": -103, "rz_map_neutron": -104}), 0),
              "normalization": (None, Real(), 1),
              "calc_freq": (None, PosReal(), 2),
              "m_file": (None, Path(), 3),
              "transform": (None, IsA(Transform, index=True), None, 0),
              "cell": (None, IsA(Cell, index=True), None)}
    prelude = (("reg", "'typ", "gap", "mgf", "trcl", "file"))
    shape = (("cell", "typ", "calc_step", "normalization", "transform", "m_file"))



class UniformElectromagneticField(PhitsObject):
    name = "uniform_electromagnetic_field"
    syntax = {"e_strength": (None, Real(), 0),
              "m_strength": (None, Real(), 1),
              "e_transform": (None, IsA(Transform, index=True), None, 0),
              "m_transform": (None, IsA(Transform, index=True), None, 0),
              "cell": (None, IsA(Cell, index=True))}
    prelude = (("reg", "elf", "mgf", "trcle", "trclm"))
    shape = (("cell", "e_strength", "m_strength", "e_transform", "m_transform"))



class MappedElectromagneticField(PhitsObject):
    name = "mapped_electromagnetic_field"
    syntax = {"typ_e": (None, FinBij({"xyz_list_charged": -1, "rz_list_charged": -2, "xyz_map_charged": -3, "rz_map_charged": -4}), 0),
              "typ_m": (None, FinBij({"xyz_list_charged": -1, "rz_list_charged": -2, "xyz_map_charged": -3, "rz_map_charged": -4}), 1),
              "calc_freq": (None, PosReal(), 2),
              "e_normalization": (None, Real(), 3),
              "m_normalization": (None, Real(), 4),
              "e_file": (None, Path(), 5),
              "m_file": (None, Path(), 6),
              "e_transform": (None, IsA(Transform, index=True), None, 0),
              "m_transform": (None, IsA(Transform, index=True), None, 0)}
    prelude = (("reg", "type", "typm", "gap", "elf", "mgf", "trcle", "trclm", "filee", "filem"))
    shape = (("cell", "typ_e", "typ_m", "gap", "e_strength", "m_strength", "e_transform", "m_transform", "e_file", "m_file"))



class DeltaRay(PhitsObject):
    name = "delta_ray"
    syntax = {"threshold": (None, PosReal(), 0),
              "cell": (None, IsA(Cell, index=True), None)}
    prelude = (("reg", "del"))
    shape = (("cell", "threshold"))



class TrackStructure(PhitsObject):
    name = "track_structure"
    syntax = {"model": (None, FinBij({None: 0, "general": -1, "optimized": 1}), 0),
              "cell": (None, IsA(Cell, index=True), None)}
    prelude = ("reg", "mID")
    shape = lambda self: (("cell", "model"))



class SuperMirror(PhitsObject):
    name = "super_mirror"
    syntax = {"reflection_surface": ((None, None), (PosReal(), PosReal()), 0),
              "material_constant": (None, Real(), 1),
              "reflectivity": (None, Real(), 2),
              "critical_q": (None, Real(), 3),
              "falloff_rate": (None, Real(), 4),
              "cutoff_width": (None, PosReal(), 5)}
    prelude = (("r-in", "r-out", "mm", "r0", "qc", "am", "wm"))
    shape = (("reflection_surface", "material_constant", "reflectivity", "critical_q", "falloff_rate", "cutoff_width"))



# refactor pending whether user Fortran routines are to be supported
# class ElasticOption(PhitsObject):
#     name = "elastic_option"
#     syntax = {""}
#     required = ["c1", "c2", "c3", "c4"]
#     positional = ["c1", "c2", "c3", "c4"]
#     optional = ["cell"]
#     shape = (("cell", "c1", "c2", "c3", "c4"))
#     prelude = (("reg", "\\c1", "\\c2", "\\c3", "\\c4"))
#     nones = {"c1": "non", "c2": "non", "c3": "non", "c4": "non"}


class FragData(PhitsObject):
    name = "frag_data"
    syntax = {"semantics": (None, FinBij({"histogram": 1, "extrapolated": 4, "interpolated": 5}), 0),
              "projectile": (None, OneOf(Particle(), Nucleide()), 1),
              "target": (None, Nucleide(), 2),
              "file": (None, Path(), 3)}
    prelude = (("opt", "proj", "targ", "'file"))
    shape = (("semantics", "projectile", "target", "file"))



class Importance(PhitsObject):
    name = "importance"
    syntax = {"particles": ("part", List(Particle()), 0),
              "importance": (None, PosReal(), 0),
              "cell": (None, IsA(Cell, index=True), None)}
    prelude = ("particles", ("reg", "imp"))
    shape = (("cell", "importance"))
    group_by = lambda self: self.particles
    separator = lambda self: self.section_title()
    max_groups = 6



class WeightWindow(PhitsObject):
    name = "weight_window"
    syntax = {"particles": ("part", List(Particle()), 0), # TODO: geometrical meshes
              "variable": (None, FinBij({"energy": "energy", "time": "time"}), 2)
              "windows": (None, List(Tuple(PosReal(), PosReal)), 1),
              "cell": (None, IsA(Cell, index=True), None)}
    prelude = lambda self: ("mesh = reg", "particles",
                            f"eng = {len(self.windows)}" if self.variable == "energy" else f"tim = {len(self.windows)}",
                            " ".join(map(lambda t: t[0], self.windows)),
                            ("reg", " ".join(f"ww{i}" for i in range(1, len(self.windows) + 1))))
    shape = lambda self: (("cell", " ".join(map(lambda t: t[1], self.windows))))
    group_by = lambda self: (self.mesh, self.particles, self.grid)
    separator = lambda self: self.section_title()
    max_groups = 6


class WWBias(PhitsObject):
    name = "ww_bias"
    syntax = {"particles": ("part", List(Particle()), 0),
              "biases": (None, List(Tuple(PosReal(), PosReal())), 1),
              "cell": (None, IsA(Cell, index=True), None)}
    prelude = lambda self: ("particles", f"eng = {len(self.biases)}", " ".join(map(lambda t: t[0], self.biases)),
                            ("reg", " ".join(f"wwb{i}" for i in range(1, len(self.biases) + 1))))
    shape = (("cell", " ".join(map(lambda t: t[1], self.biases))))
    group_by = lambda self: (self.particles, self.mesh)
    separator = lambda self: self.section_title()
    max_groups = 6




class ForcedCollisions(PhitsObject):
    name = "forced_collisions"
    syntax = {"particles": ("part", List(Particle()), 0),
              "factor": (None, PosReal(), 1),
              "force_secondaries": (None, FinBij({True: 1, False: -1}), None),
              "cell": (None, IsA(Cell, index=True), None)}

    prelude = ("particles", ("reg", "fcl"))
    shape = lambda self: (("cell", f"{self.force_secondaries * self.factor}"))

    group_by = lambda self: self.particles
    separator = lambda self: self.section_title()
    max_groups = 6



class RepeatedCollisions(PhitsObject):
    name = "repeated_collisions"
    syntax = {"particles": ("part", List(Particle()), 0),
              "collision_reps": (None, PosInt(), 1),
              "evaporation_reps":  (None, PosInt(), 2),
              "ebounds": (("emin", "emax"), (PosReal(), PosReal()), None),
              "mother": (None, List(Nulceides()), None),
              "cell": (None, IsA(Cell, index=True), None)}

    prelude = lambda self: ("particles",
                            f"mother = {len(self.mother)}" if self.mother else "",
                            " ".join(self.mother) if self.mother else "",
                            "ebounds", ("reg", "n-coll", "n-evaps"))
    shape = (("cell", "collision_reps", "evaporation_reps"))

    group_by = lambda self: (self.particles, self.mother)
    separator = lambda self: self.section_title()
    max_groups = 6




# TODO: finish this
class Multiplier(PhitsObject):
    name = "multiplier"
    syntax = {"particles": ("part", List(Particle()), 0),
              "semantics": ("interpolation", FinBij({"linear": "lin", "log": "log", "left_histogram": "glow",
                                                     "right_histogram": "ghigh"}), 1),
              "bins": (None, List(Tuple(PosReal(), PosReal())), 2)}
    shape = lambda self: (f"number = -{200 + self.index}", "interpolation", "particles", f"ne = {len(self.bins)}",
                          "\n".join(map(lambda t: f"{t[0]} {t[1]}", self.bins)))



class RegionName(PhitsObject):
    name = "region_name"
    syntax = {"reg_name": (None, Text(), 0),
              "size": (None, PosReal(), 1),
              "cell": (None, IsA(Cell, index=True), None)}
    shape = (("cell", "reg_name", "size"))



class Counter(PhitsObject):
    name = "counter"
    syntax = {"particles": ("part", List(OneOf(Particle(), Nucleide())), 0),
              "entry": (None, Integer(), None, "non"),
              "exit": (None, Integer(), None, "non"),
              "collision": (None, Integer(), None, "non"),
              "reflection": (None, Integer(), None, "non"),
              "cell": (None, IsA(Cell, index=True), None)}

    prelude = ("particles", ("reg", "in", "out", "coll", "ref"))
    shape = (("cell", "entry", "exit", "collision", "reflection"))

    group_by = lambda self: self.particles
    separator = lambda self: f"counter = {self.index}"
    max_groups = 3

class Timer(PhitsObject):
    name = "timer"
    syntax = {"entry": (None, Integer(), None),
              "exit": (None, Integer(), None),
              "collision": (None, Integer(), None),
              "reflection": (None, Integer(), None),
              "cell": (None, IsA(Cell, index=True), None)}
    prelude = (("reg", "in", "out", "coll", "ref"))
    shape = (("cell", "entry", "exit", "collision", "reflection"))
