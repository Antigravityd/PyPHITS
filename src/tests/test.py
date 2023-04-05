import sys
import os
sys.path.append(os.getcwd())

print(sys.path)

import unittest
from hypothesis import given
from hypothesis.strategies import *

from run_phits import run_phits
from base import PhitsObject
from cell import Void, Cell
from source import Cylindrical
from surface import Sphere


def test_a(cls):
    print(cls)

    req = list(map(lambda tup: tuples(*[i.strat for i in tup[1]]) if isinstance(tup[1], tuple) else tup[1].strat,
                   sorted([v for k, v in cls.syntax.items() if v[2] is not None], key=lambda tup: tup[2])))
    opt = dict(map(lambda tup: (tup[0], tuples(*[i.strat for i in tup[1][1]]) if isinstance(tup[1][1], tuple) else tup[1][1].strat),
                   [(k, v) for k, v in cls.syntax.items() if v[2] is None]))

    breakpoint()
    @given(builds(cls, *req, **opt))
    def definition_syntax_correct(ins):
        test_source = Cylindrical()
        test_surf = Sphere(1, (0, 0, 0))
        test_cell = Void(test_surf)

        if ins.name == "source":
            run_phits([ins], [test_cell], [], control="input_echo")
        elif ins.name == "cell":
            run_phits([test_source], [ins], [], control="input_echo")
        elif ins.name == "surface":
            run_phits([test_source], [Void(ins)], [], control="input_echo")
        elif ins.name == "material":
            run_phits([test_source], Cell(test_surf, material=ins, density=1.0), control="input_echo")
        elif ins.name in ["t-cross", "t-product", "t-time"]:
            run_phits([test_source], [test_cell], [ins], control="input_echo")
        elif ins.name in Cell.subobjects:
            setattr(test_cell, ins.name, ins)
            run_phits([test_source], [test_cell], [], control="input_echo")
        else:
            inp = ins.section_title() + (ins.prelude_str() if hasattr(ins, prelude) else "") + ins.definition()
            run_phits([test_source], [test_cell], [], control="input_echo", raw=inp)

    definition_syntax_correct()


    @given(builds(cls, *req, **opt), builds(cls, *req, **opt))
    def definition_bijective_wrt_equivalence(ins1, ins2): # this is Noether's first isomorphism theorem
        if ins1 == ins2:
            assert ins1.definition() == ins2.definition(), "Objects that are __eq__ should have the same definition."

        if ins1.definition() == ins2.definition():
            assert ins1 == ins2, "Objects that have the same definition should be __eq__."

    definition_bijective_wrt_equivalence()





if __name__ == '__main__':
    import base,  source, parameters #, cell, surface, tally, transform, material,
    for mod in [base, source, parameters]: # cell, misc, source, surface, tally, transform, material,  ]:
        for name, cls in mod.__dict__.items():
            if isinstance(cls, type) and cls.__module__ == mod.__name__ and issubclass(cls, PhitsObject) and cls != PhitsObject:
                test_a(cls)
