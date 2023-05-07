import sys
import os
sys.path.append(os.getcwd())

print(sys.path)

import unittest
from hypothesis import given, settings, assume, HealthCheck
from hypothesis.strategies import *

from run_phits import run_phits
from base import PhitsObject
from cell import Cell
from source import Cylindrical
from surface import Sphere
from material import Material


def test_a(cls):
    print(cls)

    req = []
    for phits_iden, valspec, idx in sorted((v for v in cls.syntax.values() if v[2] is not None), key=lambda t: t[2]):
        if isinstance(valspec, tuple):
            req.append(tuples(*[i.strat for i in valspec]))
        else:
            req.append(valspec.strat)

    opt = dict()
    for py_iden, (phits_iden, valspec, idx) in filter(lambda t: t[1][2] is None, cls.syntax.items()):
        if isinstance(valspec, tuple):
            opt[py_iden] = one_of(none(), tuples(*[i.strat for i in valspec]))
        else:
            opt[py_iden] = one_of(none(), valspec.strat)

    @composite
    def builds_right(draw, cl, re, op):
        try:
            ob = draw(builds(cl, *re, **op))
        except ValueError:
            ob = None

        assume(ob)
        return ob

    @given(builds_right(cls, req, opt))
    @settings(deadline=None, suppress_health_check=[HealthCheck.too_slow, HealthCheck.large_base_example])
    def definition_syntax_correct(ins):
        test_source = Cylindrical()
        test_surf = Sphere(1, (0, 0, 0))
        test_mat = Material([("H", 1)])
        test_cell = Cell([test_surf], test_mat, 0.5)

        if ins.name == "source":
            run_phits([test_cell], ins, [], control="output_echo_only")
        elif ins.name == "cell":
            run_phits([ins], test_source, [], control="output_echo_only")
        elif ins.name == "surface":
            run_phits([Void(ins)], test_source, [], control="output_echo_only")
        elif ins.name == "material":
            run_phits(Cell(test_surf, material=ins, density=1.0), test_source, control="output_echo_only")
        elif ins.name in ["t-cross", "t-product", "t-time"]:
            run_phits([test_cell], test_source, [ins], control="output_echo_only")
        elif ins.name in Cell.subobjects:
            setattr(test_cell, ins.name, ins)
            run_phits([test_cell], test_source, [], control="output_echo_only")
        else:
            inp = ins.section_title() + (ins.prelude_str() if hasattr(ins, "prelude") else "") + ins.definition()
            run_phits([test_cell], test_source, [], control="output_echo_only", raw=inp)

    definition_syntax_correct()


    @given(builds_right(cls, req, opt), builds_right(cls, req, opt))
    @settings(deadline=None, suppress_health_check=[HealthCheck.too_slow])
    def eq_is_definitional_equality(ins1, ins2): # this is Noether's first isomorphism theorem
        if ins1 == ins2:
            assert ins1.definition() == ins2.definition(), "Objects that are __eq__ should have the same definition."

        if ins1.definition() == ins2.definition():
            assert ins1 == ins2, "Objects that have the same definition should be __eq__."

    eq_is_definitional_equality()




if __name__ == '__main__':
    import base,  source, parameters, cell, surface, tally, transform, material, misc
    for mod in [base, source, parameters, cell, misc, source, surface, tally, transform, material]: # cell, misc, source, surface, tally, transform, material,  ]:
        for name, cls in mod.__dict__.items():
            if isinstance(cls, type) and cls.__module__ == mod.__name__ and issubclass(cls, PhitsObject) and cls != PhitsObject:
                test_a(cls)
