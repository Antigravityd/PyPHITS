from base import *
from cell import *
import collections as col


common = {"cells": ("reg", List(IsA(Cell, index=True)), 0),

          "particles": ("part", List(Particle()), None),
          }
# counters? std cutoff? bounds?
class DumpFluence(PhitsObject):
    name = "t-cross"
    syntax = {"out": (None, IsA(Cell, index=True), 0),
              "into": (None, IsA(Cell, index=True), 1),
              "area": (None, PosReal(), 2),
              }
    required = ["out", "into", "area", "data", "output_type"]
    positional = ["out", "into", "area", "data", "output_type"]
    optional = ["particles", "factor", "energy_bounds", "angle_bounds", "time_bounds"]
    shape = ((lambda self: f"{self.out.index}", lambda self: f"{self.into.index}", "area"),)
    prelude = ("mesh = reg", "particles", "unit = 1", "axis = reg", "file = cross.dmp", "factor", "output_type",
               lambda self: f"reg = {self.group_size}", lambda self: f"dump = -{len(self.data)}",
               lambda self: " ".join([str(i) for i in self.data]), ("r-out", "r-in", "'area"))
    ident_map = {"particles": "part"}
    group_by = lambda self: (self.particles, self.data, self.output_type, self.factor, self.energy_bounds, self.angle_bounds,
                             self.time_bounds)
    separator = lambda self: self.section_title()



class DumpProduct(PhitsObject):
    name = "t-cross"
    required = ["out", "into", "area", "data", "output_type"]
    positional = ["out", "into", "area", "data", "output_type"]
    optional = ["particles", "factor", "energy_bounds", "angle_bounds", "time_bounds"]
    shape = ((lambda self: f"{self.out.index}", lambda self: f"{self.into.index}", "area"),
             lambda self: f"dump = -{len(self.data)}", lambda self: " ".join([str(i) for i in self.data]))
    prelude = ("mesh = reg", "particles", "unit = 1", "axis = reg", "file = cross.dmp", "factor", "output_type",
               lambda self: f"reg = {self.group_size}")
    ident_map = {"particles": "part"}
    group_by = lambda self: (self.particles, self.data, self.output_type, self.factor, self.energy_bounds, self.angle_bounds,
                             self.time_bounds)
    separator = lambda self: self.section_title()

class DumpTime(PhitsObject):
    name = "t-cross"
    required = ["out", "into", "area", "data", "output_type"]
    positional = ["out", "into", "area", "data", "output_type"]
    optional = ["particles", "factor", "energy_bounds", "angle_bounds", "time_bounds"]
    shape = ((lambda self: f"{self.out.index}", lambda self: f"{self.into.index}", "area"),
             lambda self: f"dump = -{len(self.data)}", lambda self: " ".join([str(i) for i in self.data]))
    prelude = ("mesh = reg", "particles", "unit = 1", "axis = reg", "file = cross.dmp", "factor", "output_type",
               lambda self: f"reg = {self.group_size}")
    ident_map = {"particles": "part"}
    group_by = lambda self: (self.particles, self.data, self.output_type, self.factor, self.energy_bounds, self.angle_bounds,
                             self.time_bounds)
    separator = lambda self: self.section_title()





# Non-Dump tallies not currently supported, but they stay because we may look into converting them to the dump tallies,
# or resurrecting the Angel parsing, to support more general use-cases.
# class TrackLengthFluence(PhitsObject):
#     name = "t-track"
#     required = ["out_mesh", "units"]
#     positional = ["out_mesh", "units"]
#     optional = ["mesh_type", "outfile", "in_meshes", "dependent_var", "particles", "angel", "sangel", "title",
#                 "output_type", "type_2d", "factor", "detailed_output", "region_show", "geometry_show",
#                 "axis_titles", "epesout", "resolution", "transform", "dump", "width"]
#     shape = ("mesh_type", "out_meshes", "in_meshes", "particles", "outfile", "units", "factor", "title",
#              "angel", "sangel", "type_2d", "gshow", "rshow", "axis_titles",  "resolution", "width",
#              "transform", "epsout")
#     ident_map = {"type_2d": "2d-type", "region_show": "rshow", "geometry_show": "gshow",
#                  "dependent_var": "axis", "outfile": "file", "particles": "part", "units": "unit",
#                  "axis_titles": ("x-txt", "y-txt", "z-txt"), "transform": "trcl", "resolution": "resol",
#                  "mesh_type": "mesh", "output_type": "output"}
#     value_map = {"1/cm^2/source": 1, "1/cm^2/MeV/source": 2, "1/cm^2/Lethargy/source": 3,
#                  "cm/source": 4, "1/cm^2/nsec/source": 11, "1/cm^2/MeV/nsec/source": 12,
#                  "1/cm^2/Lethargy/nsec/source": 13, "cm/nsec/source": 14}


# class SurfaceFluence(PhitsObject):
#     name = "t-cross"
#     required = ["out_mesh", "units"]
#     positional = ["out_mesh", "units"]
#     optional = ["mesh_type", "outfile", "in_meshes", "dependent_var", "particles", "angel", "sangel", "title",
#                 "output_type", "type_2d", "factor", "detailed_output", "region_show", "geometry_show",
#                 "axis_titles", "epesout", "resolution", "transform", "dump", "width"]
#     shape = ("mesh_type", "out_meshes", "in_meshes", "particles", "outfile", "units", "factor", "title",
#              "angel", "sangel", "type_2d", "gshow", "rshow", "axis_titles",  "resolution", "width",
#              "transform", "epsout")
#     ident_map = {"type_2d": "2d-type", "region_show": "rshow", "geometry_show": "gshow",
#                  "dependent_var": "axis", "outfile": "file", "particles": "part", "units": "unit",
#                  "axis_titles": ("x-txt", "y-txt", "z-txt"), "transform": "trcl", "resolution": "resol",
#                  "mesh_type": "mesh", "output_type": "output"}
#     value_map = {"1/cm^2/source": 1, "1/cm^2/MeV/source": 2, "1/cm^2/Lethargy/source": 3,
#                  "1/cm^2/sr/source": 4, "1/cm/MeV/sr/source": 5, "1/cm^2/Lethargy/sr/source": 6,
#                  "1/cm^2/nsec/source": 11, "1/cm^2/MeV/nsec/source": 12, "1/cm^2/Lethargy/nsec/source": 13,
#                  "1/cm^2/sr/nsec/source": 14, "1/cm^2/MeV/sr/nsec/source": 15,
#                  "1/cm^2/Lethargy/sr/nsec/source": 16, "normal": 0, "x": 1, "y": 2, "z": 3}


# # class PointFluence(PhitsObject): # TODO: finish
# #     name = "t-point"
# #     required = ["out_mesh", "units"]
# #     positional = ["out_mesh", "units"]
# #     optional = ["mesh_type", "outfile", "in_meshes", "dependent_var", "particles", "angel", "sangel", "title",
# #                 "output_type", "type_2d", "factor", "detailed_output", "region_show", "geometry_show",
# #                 "axis_titles", "epesout", "resolution", "transform", "dump", "width"]
# #     shape = ("mesh_type", "out_meshes", "in_meshes", "particles", "outfile", "units", "factor", "title",
# #              "angel", "sangel", "type_2d", "gshow", "rshow", "axis_titles",  "resolution", "width",
# #              "transform", "epsout")
# #     ident_map = {"type_2d": "2d-type", "region_show": "rshow", "geometry_show": "gshow",
# #                  "dependent_var": "axis", "outfile": "file", "particles": "part", "units": "unit",
# #                  "axis_titles": ("x-txt", "y-txt", "z-txt"), "transform": "trcl", "resolution": "resol",
# #                  "mesh_type": "mesh", "output_type": "output"}
# #     value_map = {"1/cm^2/source": 1, "1/cm^2/MeV/source": 2, "1/cm^2/Lethargy/source": 3,
# #                  "1/cm^2/sr/source": 4, "1/cm/MeV/sr/source": 5, "1/cm^2/Lethargy/sr/source": 6,
# #                  "1/cm^2/nsec/source": 11, "1/cm^2/MeV/nsec/source": 12, "1/cm^2/Lethargy/nsec/source": 13,
# #                  "1/cm^2/sr/nsec/source": 14, "1/cm^2/MeV/sr/nsec/source": 15,
# #                  "1/cm^2/Lethargy/sr/nsec/source": 16, "normal": 0, "x": 1, "y": 2, "z": 3}


# class Deposition(PhitsObject):
#     name = "t-deposit"
#     required = ["out_mesh", "output_type", "units", "dependent_var"]
#     positional = ["out_mesh", "output_type", "units", "dependent_var"]
#     optional = ["mesh_type", "outfile", "in_meshes", "particles", "material", "let_material",
#                 "angel", "sangel", "title", "deposit", "type_2d",
#                 "factor", "detailed_output", "region_show", "geometry_show","axis_titles", "epsout",
#                 "resolution", "transform", "dump", "width", "gshow", "rshow", "cell"]
#     shape = (lambda self: "mesh = reg" if self.mesh_type is None else f"mesh = {self.mesh_type}",
#              lambda self: f"reg = {self.cell.index}" if self.mesh_type is None else "",
#              ("out_mesh",), ("in_meshes",), "particles", "output_type", "dependent_var",
#              lambda self: f"material = {len(self.material)}" if self.material is not None else "",
#              "material", "let_material", "outfile", "units", "factor", "title",
#              "angel", "sangel", "type_2d", "gshow", "rshow", "axis_titles",  "resolution", "width",
#              "transform", "epsout")
#     ident_map = {"type_2d": "2d-type", "region_show": "rshow", "geometry_show": "gshow",
#                  "dependent_var": "axis", "outfile": "file", "particles": "part", "units": "unit",
#                  "axis_titles": ("x-txt", "y-txt", "z-txt"), "transform": "trcl", "resolution": "resol",
#                  "mesh_type": "mesh", "output_type": "output", "let_material": "letmat"}
#     value_map = {"Gy/source": 0, "MeV/cm^3/source": 1, "MeV/source": 2, "1/source": 3,
#                  "1/nsec/source": 4, "J/m^3/source": 5, "total": 0, "per-particle": 1}
#     nones = {"outfile": "deposit.out"}



# class Time(PhitsObject):
#     name = "t-time"
# # class DepositionCorrelation(PhitsObject):
# #     name = "t-deposit2"
# #     required = ["out_mesh1", "out_mesh2", "units", "dependent_var"]
# #     positional = ["out_mesh", "output_type", "units", "dependent_var"]
# #     optional = ["mesh_type", "outfile", "in_meshes", "particles", "material", "let_material1",
# #                 "let_material2", "angel", "sangel", "title", "type_2d", "stdev_cutoff",
# #                 "factor", "axis_titles", "epsout", "volume", "cell", "other_cell"]
# #     shape = ("mesh = reg", "cell\\", ("other_cell",), "particles", "let_material1", "let_material2",
# #              "out_mesh1", "out_mesh2", "in_meshes", "units", "dependent_var", "outfile", "factor", "title",
# #              "angel", "sangel", "type_2d", "axis_titles", "epsout", "volume", "stdev_cutoff")
# #     ident_map = {"type_2d": "2d-type", "cell": "reg", "dependent_var": "axis", "outfile": "file",
# #                  "particles": "part", "units": "unit", "axis_titles": ("x-txt", "y-txt", "z-txt"),
# #                  "let_material1": "letmat1", "let_material2": "letmat2", "stdev_cutoff": "stdcut"}
# #     value_map = {"1/source": 1, "1/nsec/source": 2}
# #     nones = {"outfile": "deposit2.out"}
