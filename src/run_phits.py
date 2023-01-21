from datetime import datetime
import subprocess as sp
import os
import sys
import shutil as sh
import networkx as nx
import itertools as it
import collections as col
import re
import numpy as np
import pandas as pd
import tempfile as tf

from base import *
from cell import *
from material import *
from misc import *
from source import *
from surface import *
from tally import *
from transform import *
from tco import *


# This is complicated enough to sketch the algorithm before writing code (shocker!)
# Given a list of cells, sources, and measurements (with some wacky options thrown in), we want to generate
# all of the "background" data necessary to build those cells automatically, e.g. generate the surfaces necessary to build the cells.
# To that end, we divide the options into "cell-like," "source-like," and "tally-like," and require their specification in the initialization of the
# object of the corresponding type. This is natural: why should the surfaces I only think about in the context of one or two cells be specified separately?
# There do exist options that don't satisfy this natural classification: transforms come to mind, as they may be used in surface, source, or
# "backend" contexts.

# The user specifies a list of sources, a list of cells, and a list of tallies
# we reconstruct the input file by first going through every "PHITSable" property
# and generating a digraph of the distinct values (up to the equivalence relation of each class's __eq__), with adjacency x → y denoting
# "y needs x for its definition." This is necessary because we don't want to duplicate definitions, but we want to retain the connection
# of parent-child; we need to define two parents with __eq__ children in terms of /one/ .inp line representing the child.
# The set of values of each type are enumerated, and the graph updated to reflect each object's enumeration, representing the order in which the objects will appear in the input file.
# We can now use the child's enumeration in the definition of the parent, which is what we're after.

# In effect, we take the graph of Python object relationships, apply an algorithm that collapses __eq__ nodes into each other while retaining adjacency, (this seems like an extant problem, I just don't know its name)
# map nodes → (idx, nodes) where idx is the result of enumerating all nodes that end up in the same PHITS section, and write it all out.

# The first algorithm is, actually, best implemented by creating an explicit graph, as anything I can think of is formally equivalent to that. However, it may be space-inefficient, as the adjacency data is represented twice.
# The naïve algorithm is suboptimal: taking every node, and comparing it to every other node to construct the equivalence class explicitly is O(n^n), and doesn't use symmetry or transitivity to reduce the number of required comparisons.

# Here's a better candidate algorithm: first, initialize an empty set (since we only need lookup and insertion we can get O(1) on average; the hash values are computed via the equivalence relation in question).
# Iterate over the nodes of the graph; if the a node isn't in the set (i.e. is not equivalent to any in the set, since the hash bins are computed via __eq__), add it.
# If a node is equivalent to one in the set, change all edges that end on the node to end instead on the one to which it is equivalent.
# This is O(#nodes + connectivity)


# sources is list of tuples of a Source value and a numeric weight, cells is a list of Cell values,
# and

def make_input(cells, sources, tallies, title=str(datetime.now()), parameters=dict(), cross_sections=[], raw=""): # [Super Mirror], [Elastic Option], [Weight Window], and [WW Bias] aren't supported due to poor documentation.
    """"""
    unique = set()


    def add_to_set(an_obj, the_set, prev=None):  # Recursively add subtypes to set if they represent an "entry" in one of the sections
        if isinstance(an_obj, col.Iterable):
            for ob in an_obj:
                add_to_set(ob, the_set)
        if isinstance(an_obj, PhitsObject):
            the_set.add(an_obj)
            for name, child in an_obj.__dict__.items():
                if isinstance(child, PhitsObject):
                    the_set.add(child)
                    if child is not prev:
                        add_to_set(child, the_set, an_obj)





    add_to_set(cells, unique)
    for cell in cells: # TODO: this ideally should "just work." Every language should have ADTs
        for sur in cell.regions:
            add_to_set(sur, unique)
    add_to_set(sources, unique)
    add_to_set(tallies, unique)




    breakpoint()
    # We now have that if any two PHITS objects A and B have a property C and D (respectively) such that C == D, C /is/ D.


    # Now we break up the nodes by PHITS type; these objects will be translated into the .inp file directly in this order.
    type_divided = {"parameters": [],
                    "source": [],
                    "material": [],
                    "surface": [],
                    "cell": [],
                    "transform": [],
                    "temperature": [],
                    "mat_time_change": [],
                    "magnetic_field": [],
                    "neutron_magnetic_field": [],
                    "mapped_magnetic_field": [],
                    "uniform_electromagnetic_field": [],
                    "mapped_electromagnetic_field": [],
                    "delta_ray": [],
                    "track_structure": [],
                    "importance": [],
                    "weight_window": [],
                    "ww_bias": [],
                    "forced_collisions": [],
                    "repeated_collisions": [],
                    "volume": [],
                    "multiplier": [],
                    "mat_name_color": [],
                    "reg_name": [],
                    "counter": [],
                    "timer": [],
                    "t-track": [],
                    "t-cross": [],
                    "t-point": [],
                    "t-adjoint": [],
                    "t-deposit": [],
                    "t-deposit2": [],
                    "t-heat": [],
                    "t-yield": [],
                    "t-product": [],
                    "t-dpa": [],
                    "t-let": [],
                    "t-sed": [],
                    "t-time": [],
                    "t-interact": [],
                    "t-dchain": [],
                    "t-wwg": [],
                    "t-wwbg": [],
                    "t-volume": [],
                    "t-gshow": [],
                    "t-rshow": [],
                    "t-3dshow": []}

    for node in unique:
        type_divided[node.name].append(node)


    for section, entries in type_divided.items():
        for idx, value in enumerate(entries):
            value.index = idx+1 # this allows us to access the position in which the value will appear if given value alone.

    representatives = {n: n for n in it.chain.from_iterable(type_divided.values())}

    def adjust_subobjects(an_obj, dic, prev=None): # Recursively replace redundant subtypes with the representative in the dict
        if isinstance(an_obj, tuple):
            for ob in an_obj:
                adjust_subobjects(ob, dic)
        elif isinstance(an_obj, PhitsObject):
            for name, child in an_obj.__dict__.items():
                if isinstance(child, PhitsObject):
                    representative = representatives[child]
                    setattr(an_obj, name, representative)
                    if child is not prev:
                        adjust_subobjects(child, dic, an_obj)

    adjust_subobjects(cells, representatives)
    for cell in cells:
        for i, sur in enumerate(cell.regions):
            adjust_subobjects(sur, representatives)
            cell.regions = tuple(representatives[s] for s in cell.regions)
    adjust_subobjects(sources, representatives)
    adjust_subobjects(tallies, representatives)




    inp = ""
    def add_defs(obj_type):
        nonlocal inp
        if obj_type in type_divided:
            if type_divided[obj_type]:
                objs = type_divided[obj_type]
                if hasattr(objs[0], "group_by"):
                    grouped = it.groupby(objs, objs[0].group_by)
                    assert len(grouped) <= objs[0].max_groups, ValueError(f"Too many {obj_type} groups.")
                    for group in grouped:
                        rep = group[0]
                        inp += rep.separator(rep)
                        if hasattr(rep, "prelude_shape"):
                            inp += rep.prelude()
                        for obj in group:
                            inp += obj.definition()
                else:
                    inp += objs[0].section_title()
                    if hasattr(objs[0], "prelude_shape"):
                        inp += objs[0].prelude()
                    for obj in objs:
                        inp += obj.definition()



    breakpoint()
    inp += "[Title]\n"
    inp += title + '\n'

    if parameters or any(param.empty() for param in type_divided["parameters"]):
        inp += "[Parameters]\n"
        for var, val in parameters.items(): # directly passed global parameters
            if var not in {"totfact", "iscorr"}:
                inp += f"{var} = {val}\n"

        add_defs("parameters") # parameters associated with object declarations, but that need to be in this global context


    inp += "[Source]\n"
    if "totfact" in parameters:
        val = parameters["totfact"]
        inp += f"totfact = {val}"
    if "iscorr" in parameters:
        val = parameters["iscorr"]
        inp += f"iscorr = {val}"

    if isinstance(sources, col.Iterable):
        for source, weight in sources:
            inp += f"<source> = {weight}\n"
            inp += source.definition()
    else:
        inp += sources.definition()


    add_defs("material")
    add_defs("surface")
    add_defs("cell")
    add_defs("transform")
    add_defs("mat_time_change")
    add_defs("magnetic_field")
    add_defs("neutron_magnetic_field")
    add_defs("mapped_magnetic_field")
    add_defs("uniform_electromagnetic_field")
    add_defs("mapped_electromagnetic_field")
    add_defs("delta_ray")
    add_defs("track_structure")

    if cross_sections:
        inp += "[Frag Data]\n"
        inp += "opt proj targ file\n"
        for dic in cross_sections:
            opt = "opt"
            proj = "proj"
            targ = "targ"
            filename = "file"
            inp += f"{dic[opt]} {dic[proj]} {dic[targ]} {dic[filename]}\n"

    if type_divided["importance"]:
        if len(type_divided["importance"]) > 6:
            raise ValueError("More than 6 [Importance] sections.")
        for imps in it.groupby(type_divided["importance"], lambda imp: imp.particles):
            inp += "[Importance]\n"
            inp += f"part = {imps[0].particles}\n"
            inp += "reg imp\n"
            for imp in imps:
                inp += f"{imp.cell.index} {imp.importance}\n"

  # if type_divided["weight_window"]:
  # if type_divided["ww_bias"]:

    if type_divided["forced_collisions"]:
        if len(type_divided["forced_collisions"]) > 6:
            raise ValueError("More than 6 [Forced Collision] sections.")
        for fcls in it.groupby(type_divided["forced_collisions"], lambda fcl: fcl.particles):
            inp += "[Forced Collisions]\n"
            inp += f"part = {fcls[0].particles}\n"
            inp += "reg imp\n"
            for fcl in fcls:
                inp += f"{fcl.cell.index} {fcl.fcl}\n"

  # if type_divided["repeated_collisions"]:
  # if type_divided["multiplier"]:

    add_defs("mat_name_color")
    add_defs("reg_name")

  # if type_divided["counter"]

    add_defs("timer")

    add_defs("t-track")
    add_defs("t-cross")
    add_defs("t-point")
    add_defs("t-deposit")
    # ...

    inp += raw

    return inp




def capture_result(return_type): # -> pandas.DataFrame | numpy.array | dict
    files = []
    for name in os.listdir():
        if re.search("\.out$", name) is not None:
            files.append(name)


    breakpoint()
    output = []
    for fil in files:
        with open(fil, "r") as f:
            schema = []
            data = dict()
            in_data = False
            for line in f:
                if line[0] == "#":
                    continue
                if in_data:
                    if line in ['\n', '\r\n']:
                        in_data = False
                    else:
                        for idx, datum in enumerate(line.split()):
                            data[idx].append(float(datum))
                else:
                    if re.search("^H.?:", line, flags=re.IGNORECASE) is not None:
                        in_data = True
                        line = re.sub("^H.?:", "", line, flags=re.IGNORECASE)
                        sep  = line.split()
                        for idx, col in enumerate(sep):
                            schema.append(col)
                            data[idx] = []

            breakpoint()
            for idx, colname in enumerate(schema):
                data[colname] = data[idx]
                del data[idx]

            if return_type == "dict":
                output.append(data)
            elif return_type == "pandas":
                output.append(pd.DataFrame(data))
            elif return_type == "numpy":
                array = []
                for col, lst in data.items():
                    array.append(lst)
                output.append(np.array(array))

    return output



class TempWorkingDir(): # this is probably all kinds of racy and unsafe, but eh
    def __init__(self, name):
        self.name = name
    def __enter__(self):
        self.directory = os.mkdir(self.name)
        os.chdir(self.name)
    def __exit__(self, a, b, c):
        os.chdir("..")
        sh.rmtree(self.name)

def run_phits(sources, cells, tallies, command="phits", throws=False, filename="phits", return_type="dict", **kwargs):
    # TODO: consider how to read stdout/output files into returnable formats
    # WARNING: setting the command variable opens up shell injection attacks, as sp.run() with
    # shell=True is done unfiltered. Should see about using shlex.quote() to sanitize, since title may be specified by the user

    with TempWorkingDir("temp_PHITS") as newdir:
        inp = make_input(sources, cells, tallies, **kwargs)
        print(inp)
        with open(f"{filename}.inp", "w") as inp_file:
            inp_file.write(inp)

        output = sp.run(f"phits '{filename}.inp'", shell=True, capture_output=True, text=True, check=throws)
        print(output)


        result = capture_result(return_type)

        return result




# def parse_input(file_handle) -> :#< tuple of all arguments of run_phits function>
