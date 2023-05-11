import sys
from base import *
from transform import *
from surface import surface_spec
import itertools as it

# no temperature; do that at the cell level for now

class MagneticField(PhitsObject): # Right now, the only way to set this is to do setattr(<cell>, "magnetic_field", <MagneticField>) since one cannot pass the cell to create the magnetic field while the cell is being initialized.
    name = "magnetic_field"
    syntax = {"typ": (None, FinBij({"dipole": 2, "quadrupole": 4}), 0),
              "strength": (None, Real(), 1),
              "calc_freq": (None, PosReal(), None, 0.0),
              "transform": (None, IsA(Transform, index=True), None, 0),
              "time": (None, PosReal(), None, "non"),
              }
    superobjects = ["cell"]
    prelude = (("reg", "'typ", "'gap", "mgf", "trcl", "'time"),)
    shape = (("cell", "typ", "calc_freq", "strength", "transform", "time"),)

    group_by = lambda self: type(self).__name__
    separator = lambda self: self.section_title()

class NeutronMagneticField(PhitsObject):
    name = "magnetic_field"
    syntax = {"typ": (None, FinBij({"identified": 60, "nograv": 61, "dipole": 102,
                                    "quadrupole": 104, "sextupole": 106}), 0),
              "strength": (None, Real(), 1),
              "calc_freq": (None, PosReal(), None, 0.0),
              "polarization": (None, Real(), None, "non"),
              "transform": (None, IsA(Transform, index=True), None, 0),
              "time": (None, PosReal(), None, "non"),
              }
    superobjects = ["cell"]
    prelude = (("reg", "'typ", "'gap", "mgf", "trcl", "polar", "'time"),)
    shape = (("cell", "typ", "calc_freq", "strength", "transform", "polarization", "time"),)

    group_by = lambda self: type(self).__name__
    separator = lambda self: self.section_title()


class MappedMagneticField(PhitsObject):
    name = "magnetic_field"
    syntax = {"typ": (None, FinBij({"xyz_list_charged": -1, "rz_list_charged": -2, "xyz_map_charged": -3, "rz_map_charged": -4,
                                    "xyz_list_neutron": -101, "rz_list_neutron": -102, "xyz_map_neutron": -103, "rz_map_neutron": -104}), 0),
              "normalization": (None, Real(), 1),
              "calc_freq": (None, PosReal(), 2),
              "m_file": (None, Path(), 3),
              "transform": (None, IsA(Transform, index=True), None, 0),
              }
    superobjects = ["cell"]
    prelude = (("reg", "'typ", "gap", "mgf", "trcl", "file"),)
    shape = (("cell", "typ", "calc_freq", "normalization", "transform", "m_file"),)

    group_by = lambda self: type(self).__name__
    separator = lambda self: self.section_title()



class UniformElectromagneticField(PhitsObject):
    name = "electromagnetic_field"
    syntax = {"e_strength": (None, Real(), 0),
              "m_strength": (None, Real(), 1),
              "e_transform": (None, IsA(Transform, index=True), None, 0),
              "m_transform": (None, IsA(Transform, index=True), None, 0),
              }
    superobjects = ["cell"]
    prelude = (("reg", "elf", "mgf", "trcle", "trclm"),)
    shape = (("cell", "e_strength", "m_strength", "e_transform", "m_transform"),)

    group_by = lambda self: type(self).__name__
    separator = lambda self: self.section_title()



class MappedElectromagneticField(PhitsObject):
    name = "electromagnetic_field"
    syntax = {"typ_e": (None, FinBij({"xyz_list_charged": -1, "rz_list_charged": -2, "xyz_map_charged": -3, "rz_map_charged": -4}), 0),
              "typ_m": (None, FinBij({"xyz_list_charged": -1, "rz_list_charged": -2, "xyz_map_charged": -3, "rz_map_charged": -4}), 1),
              "calc_freq": (None, PosReal(), 2),
              "e_normalization": (None, Real(), 3),
              "m_normalization": (None, Real(), 4),
              "e_file": (None, Path(), 5),
              "m_file": (None, Path(), 6),
              "e_transform": (None, IsA(Transform, index=True), None, 0),
              "m_transform": (None, IsA(Transform, index=True), None, 0)}
    superobjects = ["cell"]
    prelude = (("reg", "type", "typm", "gap", "elf", "mgf", "trcle", "trclm", "filee", "filem"),)
    shape = (("cell", "typ_e", "typ_m", "calc_freq", "e_normalization", "m_normalization", "e_transform", "m_transform", "e_file", "m_file"),)

    group_by = lambda self: type(self).__name__
    separator = lambda self: self.section_title()



class DeltaRay(PhitsObject):
    name = "delta_ray"
    syntax = {"threshold": (None, RealBetween(1, None), 0),
              }
    superobjects = ["cell"]
    prelude = (("reg", "del"),)
    shape = (("cell", "threshold"),)



class TrackStructure(PhitsObject):
    name = "track_structure"
    syntax = {"model": (None, FinBij({"none": 0, "general": -1, "optimized": 1}), 0)}
    superobjects = ["cell"]
    prelude = (("reg", "mID"),)
    shape = lambda self: (("cell", "model"),)



# TODO: seemingly irremovable circularity here (it wants cells, not surfaces for reflection_surface)
# class SuperMirror(PhitsObject):
#     name = "super_mirror"
#     syntax = {"reflection_surface": ((None, None), (surface_spec, surface_spec), 0),
#               "material_constant": (None, Real(), 1),
#               "reflectivity": (None, Real(), 2),
#               "critical_q": (None, Real(), 3),
#               "falloff_rate": (None, Real(), 4),
#               "cutoff_width": (None, PosReal(), 5)}
#     superobjects = ["cell"]
#     prelude = (("r-in", "r-out", "mm", "r0", "qc", "am", "wm"),)
#     shape = (("reflection_surface", "material_constant", "reflectivity", "critical_q", "falloff_rate", "cutoff_width"),)



class ElasticOption(PhitsObject):
    name = "elastic_option"
    syntax = {"c1": (None, PosReal(), 1),
              "c2": (None, PosReal(), 2),
              "c3": (None, PosReal(), 3),
              "c4": (None, PosReal(), 4)}

    prelude = (("reg", "'c1", "'c2", "'c3", "'c4"),)
    shape = (("cell", "c1", "c2", "c3", "c4"),)
    superobjects = ["cell"]



class FragData(PhitsObject):
    name = "frag_data"
    syntax = {"semantics": (None, FinBij({"histogram": 1, "extrapolated": 4, "interpolated": 5}), 0),
              "projectile": (None, OneOf(FinBij({"proton": "proton", "neutron": "neutron"}), Nuclide(fake=True)), 1),
              "target": (None, Nuclide(fake=True), 2),
              "file": (None, Path(), 3)}
    superobjects = ["cell"]
    prelude = (("opt", "proj", "targ", "'file"),)
    shape = (("semantics", "projectile", "target", "file"),)



class Importance(PhitsObject):
    name = "importance"
    syntax = {"particles": ("part", List(Particle(), unique=True), 0),
              "importance": (None, PosReal(), 1),
              }
    superobjects = ["cell"]
    prelude = ("particles", ("reg", "imp"))
    shape = (("cell", "importance"),)
    group_by = lambda self: self.particles
    separator = lambda self: self.section_title()
    max_groups = 6

    @classmethod
    def global_restrictions(self, type_divided):
        all_particles = list(it.chain.from_iterable(map(lambda x: x.particles, type_divided["importance"])))
        if len(set(all_particles)) < len(all_particles):
            raise ValueError("Integration problem: all Importances must have mutually disjoint lists of particles.")


class WeightWindow(PhitsObject):
    name = "weight_window"
    syntax = {"particles": ("part", List(Particle(), unique=True), 0),
              "variable": (None, FinBij({"energy": "energy", "time": "time"}), 2),
              "windows": (None, List(Tuple(PosReal(), PosReal())), 1),
              }
    superobjects = ["cell"]
    prelude = lambda self: ("mesh = reg", "particles",
                            f"eng = {len(self.windows)}" if self.variable == "energy" else f"tim = {len(self.windows)}",
                            " ".join(map(lambda t: str(t[0]), self.windows)),
                            ("reg", " ".join(f"ww{i}" for i in range(1, len(self.windows) + 1))))
    shape = lambda self: (("cell", " ".join(map(lambda t: str(t[1]), self.windows))),)
    group_by = lambda self: (self.particles, self.variable)
    separator = lambda self: self.section_title()
    max_groups = 6

    @classmethod
    def global_restrictions(self, type_divided):
        all_particles = list(it.chain.from_iterable(map(lambda x: x.particles, type_divided["weight_window"])))
        if len(set(all_particles)) < len(all_particles):
            raise ValueError("Integration problem: all WeightWindows must have mutually disjoint lists of particles.")



class WWBias(PhitsObject):
    name = "ww_bias"
    syntax = {"particles": ("part", List(Particle(), unique=True), 0),
              "biases": (None, List(Tuple(PosReal(), PosReal())), 1),
              }
    superobjects = ["cell"]
    prelude = lambda self: ("particles", f"eng = {len(self.biases)}", " ".join(map(lambda t: str(t[0]), self.biases)),
                            ("reg", " ".join(f"wwb{i}" for i in range(1, len(self.biases) + 1))))
    shape = lambda self: (("cell", " ".join(map(lambda t: str(t[1]), self.biases))),)
    group_by = lambda self: self.particles
    separator = lambda self: self.section_title()
    max_groups = 6




class ForcedCollisions(PhitsObject):
    name = "forced_collisions"
    syntax = {"particles": ("part", List(Particle(), unique=True), 0),
              "factor": (None, RealBetween(-1, 1), 1),
              "force_secondaries": (None, FinBij({True: 1, False: -1}), None),
              }

    superobjects = ["cell"]
    prelude = ("particles", ("reg", "fcl"))
    shape = lambda self: (("cell", f"{self.force_secondaries * self.factor}" if self.force_secondaries is not None \
                           else str(self.factor)),)

    group_by = lambda self: self.particles

    separator = lambda self: self.section_title()
    max_groups = 6

    def restrictions(self):
        if "electron" in self.particles or "positron" in self.particles:
            raise ValueError(f"ForcedCollision does not accept electrons or positrons; got {self.particles}.")

    def global_restrictions(self, type_divided):
        all_particles = list(it.chain.from_iterable(map(lambda x: x.particles, type_divided["forced_collisions"])))
        if len(set(all_particles)) < len(all_particles):
            raise ValueError("Integration problem: all ForcedCollisions must have mutually disjoint lists of particles.")



class RepeatedCollisions(PhitsObject):
    name = "repeated_collisions"
    syntax = {"particles": ("part", List(Particle(fake=True), unique=True), 0),
              "collision_reps": (None, PosInt(), 1),
              "evaporation_reps":  (None, PosInt(), 2),
              "mother": (None, List(Nuclide(fake=True)), 3),
              "ebounds": (("emin", "emax"), (PosReal(), PosReal()), None),

              }

    superobjects = ["cell"]
    prelude = lambda self: ("particles",
                            f"mother = {len(self.mother)}" if self.mother else "",
                            ("mother",),
                            "ebounds", ("reg", "n-coll", "n-evap"))
    shape = (("cell", "collision_reps", "evaporation_reps"),)

    group_by = lambda self: (self.particles, self.mother)
    separator = lambda self: self.section_title()
    max_groups = 6

    def restrictions(self):
        if self.collision_reps * self.evaporation_reps <= 1 or self.collision_reps * self.evaporation_reps >= 2_147_483_647:
            raise ValueError(f"RepeatedCollisions' product of repititions must be more than 1 as an int32;"
                             f" got collsion_reps={self.collision_reps} and evaporation_reps={self.evaporation_reps}.")
        if self.ebounds is not None and self.ebounds[0] >= self.ebounds[1]:

            raise ValueError(f"RepeatedCollisions' ebounds must be ordered; got {self.ebounds}.")

    @classmethod
    def global_restrictions(self, type_divided):
        # Doesn't work because cell isn't set at __init__
        for rc in type_divided["repeated_collisions"]:
            possible = set(map(lambda x: kf_encode(x[0]), rc.cell.material.composition))
            if any(kf_encode(x) not in possible for x in rc.mother):
                raise ValueError(f"Integration problem: RepeatedCollisions' mother nuclei must be among its cell's material's nuclei;"
                                 f" got {set(rc.mother) - possible} extra.")



class Multiplier(PhitsObject):
    name = "multiplier"
    syntax = {"particles": ("part", List(Particle(), unique=True, max_len=6), 0),
              "semantics": ("interpolation", FinBij({"linear": "lin", "log": "log", "left_histogram": "glow",
                                                     "right_histogram": "ghigh"}), 1),
              "bins": (None, List(Tuple(PosReal(), PosReal())), 2)}
    shape = lambda self: (f"number = -{200 + self.index}", "semantics", "particles", f"ne = {len(self.bins)}",
                          "\n".join(map(lambda t: f"{t[0]} {t[1]}", self.bins)))



class RegionName(PhitsObject):
    name = "reg_name"
    syntax = {"reg_name": (None, Text(), 0),
              "size": (None, PosReal(), 1),
              }
    superobjects = ["cell"]
    shape = (("cell", "reg_name", "size"),)



# TODO: optional arguments?
class Counter(PhitsObject):
    name = "counter"
    syntax = {"particles": ("part", List(OneOf(Particle(), Nuclide()), max_len=20, unique=True), 0),
              "entry": (None, Between(-9999, 10000), 1),
              "exit": (None, Between(-9999, 10000), 2),
              "collision": (None, Between(-9999, 10000), 3),
              "reflection": (None, Between(-9999, 10000), 4),
              }

    superobjects = ["cell"]
    prelude = ("particles", ("reg", "in", "out", "coll", "ref"))
    shape = (("cell", "entry", "exit", "collision", "reflection"),)

    group_by = lambda self: self.particles
    separator = lambda self: self.section_title() + f"counter = {self.index}\n"
    max_groups = 3



class Timer(PhitsObject):
    name = "timer"
    syntax = {"entry": (None, FinBij({"zero": -1, "nothing": 0, "stop": 1}), 1),
              "exit": (None, FinBij({"zero": -1, "nothing": 0, "stop": 1}), 2),
              "collision": (None, FinBij({"zero": -1, "nothing": 0, "stop": 1}), 3),
              "reflection": (None, FinBij({"zero": -1, "nothing": 0, "stop": 1}), 4),
              }
    superobjects = ["cell"]
    prelude = (("reg", "in", "out", "coll", "ref"),)
    shape = (("cell", "entry", "exit", "collision", "reflection"),)

__pdoc__ = dict()
__pdoc__["builds"] = False
__pdoc__["slices"] = False
for name, cl in list(sys.modules[__name__].__dict__.items()):
    if type(cl) == type and issubclass(cl, PhitsObject) and cl != PhitsObject:
        __pdoc__[cl.__name__] = cl.__doc__ + cl.syntax_desc() if cl.__doc__ else cl.syntax_desc()
