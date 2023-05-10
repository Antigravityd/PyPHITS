from datetime import datetime
import subprocess as sp
import itertools as it
from functools import reduce
import collections as col
import numpy as np
import pandas as pd
import tempfile as tf
import re
import os

from base import PhitsObject
from parameters import Parameters
from cell import OuterVoid
from dmp_reader import read_dump


def _max_line_len(inp: str) -> int:
    return max(map(len, inp.split("\n")))


def make_input(cells, sources, tallies, title: str = str(datetime.now()), cross_sections=[], multipliers=[],
               raw="", outer_void_properties: dict = dict(), **kwargs) -> str:
    """Given a situation, produces a corresponding input file.

    Required arguments:

    | Name | Position | Description |
    | ---- | -------- | ----------- |
    | cells | 0 | A list of `PhitsObject`s with `name == "cell"`.|
    | sources | 1 | Either a single `PhitsObject` with `name == "source"`, or a list of tuples (<source object>, <weight>).|
    | tallies | 2 | A list of objects of type `DumpFluence`, `DumpProduction`, or `DumpTime`.|

    Optional arguments:

    | Name | Description |
    | ---- | ----------- |
    | title | A string to paste in the `[Title]` section. |
    | parameters | Some globally-passed options, fed directly into a `Parameters` object. |
    | cross_sections | A list of `FragData` objects. |
    | raw | A string that's appended to the end of the .inp---do unsupported stuff manually with this.|
    | kwargs | Anything extra is used to create a Parameters() object. |
    """

    # Problem: you can have different objects that are "essentially the same" appearing in the object tree.
    # Solution: exploit the __eq__ and __hash__ defined on PhitsObject
    unique = set()

    def add_to_set(an_obj, the_set, prev=None):  # Recursively add subtypes to set if they represent an "entry" in one of the sections
        if isinstance(an_obj, col.Iterable):
            for ob in an_obj:
                if isinstance(ob, PhitsObject) and ob is not prev:
                    add_to_set(ob, the_set)
        if isinstance(an_obj, PhitsObject):
            the_set.add(an_obj)
            for name, child in an_obj.__dict__.items():
                if child is not prev:
                    add_to_set(child, the_set, an_obj)

    add_to_set(cells, unique)
    add_to_set(sources, unique)
    add_to_set(tallies, unique)
    add_to_set(cross_sections, unique)
    add_to_set(multipliers, unique)



    # We now have that if any two PHITS objects A and B have attributes C and D (respectively) such that C == D, C /is/ D.


    # Problem: before this function is invoked, we can't give objects an ID number by which they're referenced in the .inp---so
    # they don't have IDs yet.
    # Solution: put all the objects in an indexed structure; index + 1 := ID.
    type_divided = {"parameters": [],
                    "source": [],
                    "material": [],
                    "surface": [],
                    "cell": [],
                    "transform": [],
                    "temperature": [],
                    "mat_time_change": [],
                    "magnetic_field": [],
                    "electromagnetic_field": [],
                    "delta_ray": [],
                    "track_structure": [],
                    # "super_mirror": [],
                    "elastic_option": [],
                    # "data_max": [],
                    "frag_data": [],
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
                    # "t-track": [],
                    "t-cross": [],
                    # "t-point": [],
                    # "t-adjoint": [],
                    # "t-deposit": [],
                    # "t-deposit2": [],
                    # "t-heat": [],
                    # "t-yield": [],
                    "t-product": [],
                    # "t-dpa": [],
                    # "t-let": [],
                    # "t-sed": [],
                    "t-time": [],
                    # "t-interact": [],
                    # "t-dchain": [],
                    # "t-wwg": [],
                    # "t-wwbg": [],
                    # "t-volume": [],
                    # "t-gshow": [],
                    # "t-rshow": [],
                    # "t-3dshow": []
                    }
    if kwargs:
        type_divided["parameters"].append(Parameters(**kwargs))
    for node in unique:
        type_divided[node.name].append(node)

    toset = OuterVoid([], **outer_void_properties)
    toset.regions = (~reduce(lambda c1, c2: c1 | c2, type_divided["cell"])).regions

    type_divided["cell"].append(toset)


    for section, entries in type_divided.items():
        for idx, value in enumerate(entries):
            value.index = idx+1

    # Problem: while we've chosen a set of representatives for equivalence classes under PhitsObject.__eq__, the objects themselves
    # don't have subobjects with index attributes pointing to the representative---there will be None showing up all over the output.
    # Solution: replace all members of an equivalence class in the object tree with their representative (whose index is defined above).
    representatives = {n: n for n in it.chain.from_iterable(type_divided.values())} # necessary because `unique` doesn't have idx

    def adjust_subobjects(an_obj, ason=(None, None)): # Recursively replace redundant subtypes with the representative in the dict
        if isinstance(an_obj, col.Iterable):
            for ob in an_obj:
                if isinstance(ob, PhitsObject) and ob is not ason[1]:
                    adjust_subobjects(ob)

        elif isinstance(an_obj, PhitsObject):
            if ason != (None, None):
                representative = representatives[an_obj]
                setattr(ason[1], ason[0], representative)
            for name, child in an_obj.__dict__.items():
                if child is not ason[1]:
                    adjust_subobjects(child, ason=(name, an_obj))

    adjust_subobjects(cells)
    adjust_subobjects(sources)
    adjust_subobjects(tallies)
    adjust_subobjects(cross_sections)
    adjust_subobjects(multipliers)


    # Now, we can make the input file.
    inp = ""
    def add_defs(obj_type):
        nonlocal inp
        if obj_type in type_divided:
            if type_divided[obj_type]:
                objs = type_divided[obj_type]
                type_rep = objs[0]
                if hasattr(type_rep, "group_by") and callable(type_rep.group_by):
                    grouped = [(k, list(v)) for k, v in it.groupby(sorted(objs, key=lambda x: x.group_by()), lambda x: x.group_by())]
                    if hasattr(type_rep, "max_groups") and type_rep.max_groups is not None:
                        assert len(grouped) <= type_rep.max_groups, ValueError(f"Too many {obj_type} groups.")
                    for key, group in grouped:
                        group = list(group)
                        inp += group[0].separator()
                        gs = len(group)
                        for obj in group:
                            obj.group_size = gs
                        if hasattr(group[0], "prelude"):
                            inp += group[0].prelude_str()
                        for obj in group:
                            inp += obj.definition()
                else:
                    inp += type_rep.section_title()
                    if hasattr(type_rep, "prelude"):
                        inp += type_rep.prelude_str()
                    for obj in objs:
                        inp += obj.definition()



    inp += "[Title]\n"
    inp += title + '\n'

    if any(not param.empty() for param in type_divided["parameters"]):
        add_defs("parameters") # parameters associated with object declarations, but that need to be in this global context.
        # for var, val in parameters.items(): # directly passed global parameters
        #     if var not in {"totfact", "iscorr"}: # TODO: document these two
        #         inp += f"{var} = {val}\n"




    inp += "[Source]\n"
    if "totfact" in kwargs:
        val = kwargs["totfact"]
        inp += f"totfact = {val}\n"
    if "iscorr" in kwargs:
        val = kwargs["iscorr"]
        inp += f"iscorr = {val}\n"

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
    add_defs("electromagnetic_field")
    add_defs("delta_ray")
    add_defs("track_structure")
    # add_defs("super_mirror")
    add_defs("elastic_option")
    # add_defs("data_max")
    add_defs("frag_data")
    add_defs("importance")
    add_defs("weight_window")
    add_defs("ww_bias")
    add_defs("forced_collisions")
    add_defs("repeated_collisions")
    add_defs("multiplier")
    add_defs("mat_name_color")
    add_defs("reg_name")
    add_defs("counter")
    add_defs("timer")
    add_defs("t-cross")
    add_defs("t-product")
    add_defs("t-time")

    inp += raw

    if _max_line_len(inp) > 200:
        raise RuntimeError("PHITS line limit reached.")
    else:
        return inp





def run_phits(cells, sources, tallies, command: str = "phits", hard_error: bool = True, filename: str = "phits.inp",
              return_type: str = "dict", **make_input_kwargs):
    """Given a scenario, calls `make_input` to generate a corresponding input file, runs PHITS on it in a temporary directory,
    and then collects and returns the resulting output as nice Python objects.

    Required arguments:

    | Name | Position | Description |
    | ---- | -------- | ----------- |
    | `cells` | 0 | A list of `PhitsObject`s with `name == "cell"`.|
    | `sources` | 1 | Either a single `PhitsObject` with `name == "source"`, or a list of tuples (<source object>, <weight>).|
    | `tallies` | 2 | A list of objects of type `DumpFluence`, `DumpProduction`, or `DumpTime`.|

    Optional arguments:

    | Name | Description |
    | ---- | ----------- |
    | `command` | The shell command to invoke on the generated file. |
    | `hard_error` | If truthy, raise an error and halt if PHITS encounters one. Otherwise, simply print the error and continue
    (helpful in avoiding "I let it run all night and it crashed the minute I left the room" scenarios). |
    | `filename` | The name of the input file on which to call PHITS. Of little utility except for debugging. |
    | `return_type` | Either "dict", "numpy", or "pandas", corresponding to the desired result format. |
    """
    with tf.TemporaryDirectory() as newdir:
        inp = make_input(cells, sources, tallies, **make_input_kwargs)
        name = os.path.join(newdir, filename)
        with open(os.path.join(newdir, filename), "w") as inp_file:
            inp_file.write(inp)

        try:
            # TODO: PHITS actualy **DOESN'T FKING SET EXIT CODES ON ERROR** so will have to grep the output for "Error"...
            out = sp.run(["phits", filename], capture_output=True, text=True, cwd=newdir)
            assert not re.search("(?i:Error)", out.stdout), "PHITS Error." # this REALLY sucks. Thank PHITS.
            assert out.returncode == 0, "PHITS Error"
        except AssertionError as error:
            r = f"PHITS exited with code {out.returncode}.\n"
            r += f"stdout: {out.stdout}\n"
            r += f"stderr: {out.stderr}\n"
            r += "Offending input file:\n"
            for idx, line in enumerate(inp.split("\n")):
                r += f"{idx}:    {line}\n"

            if hard_error:
                raise RuntimeError(r)
            else:
                print(r)

        result = dict()
        for t in tallies:

            dfile = ""
            if t.name == "t-cross":
                dfile = os.path.join(newdir, f"cross{t.index}_dmp")
            elif t.name == "t-product":
                dfile = os.path.join(newdir, f"product{t.index}_dmp")
            elif t.name == "t-time":
                dfile = os.path.join(newdir, f"time{t.index}_dmp")
            breakpoint()
            result[t] = read_dump(dfile, t.data, return_type)


        return result
