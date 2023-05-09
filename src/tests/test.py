import sys
import os
import re
sys.path.append(os.getcwd())



from hypothesis import given, settings, assume, HealthCheck, Phase
from hypothesis.strategies import *


from run_phits import run_phits, make_input
from base import PhitsObject
from cell import Cell, Void
from source import Cylindrical
from surface import Sphere
from material import Material
from distribution import EnergyDistribution


def test_a(cls):
    print(cls)


    req = []
    for phits_iden, valspec, idx in sorted((v for v in cls.syntax.values() if v[2] is not None), key=lambda t: t[2]):
        if isinstance(valspec, tuple):
            req.append(tuples(*[i.strat for i in valspec]))
        else:
            req.append(valspec.strat)

    opt = dict()
    for py_iden, (phits_iden, valspec, idx, *s) in filter(lambda t: t[1][2] is None, cls.syntax.items()):
        if isinstance(valspec, tuple):
            opt[py_iden] = one_of(none(), tuples(*[i.strat for i in valspec]))
        else:
            opt[py_iden] = one_of(none(), valspec.strat)

    # Some things are more troublesome to generate than there worth, e.g. all of the acceptable element-isotope combinations
    if cls.name == 'source':
        if "mask" in opt:
            del opt["mask"]



    @composite
    def builds_right(draw, cl, re, op):
        try:
            ob = draw(builds(cl, *re, **op))
        except ValueError:
            ob = None

        assume(ob)
        return ob

    @given(builds_right(cls, req, opt))
    @settings(deadline=None, suppress_health_check=[HealthCheck.too_slow, HealthCheck.large_base_example,
                                                    HealthCheck.data_too_large],
              phases=(Phase.explicit, Phase.reuse, Phase.generate, Phase.target,
                      # Phase.shrink
                      ))
    def definition_syntax_correct(ins):
        test_source = Cylindrical(["1H"], EnergyDistribution(bins=[(1.0, 1.0, 1.0)]))
        test_surf = Sphere(1, (0, 0, 0))
        test_mat = Material([("H", 1)])
        test_cell = Cell([test_surf], test_mat, 0.5)

        try:
            if ins.name == "source":
                inp = make_input([test_cell], ins, [], control="output_echo_only")
                assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                run_phits([test_cell], ins, [], control="output_echo_only")
            elif ins.name == "cell":
                inp = make_input([ins, test_cell], test_source, [], control="output_echo_only")
                assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                run_phits([ins, test_cell], test_source, [], control="output_echo_only")

            elif ins.name == "surface":
                inp = make_input([Void([ins]), test_cell], test_source, [], control="output_echo_only")
                assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                run_phits([Void([ins]), test_cell], test_source, [], control="output_echo_only")
            elif ins.name == "material":
                inp = make_input(Cell(test_surf, material=ins, density=1.0), test_source, control="output_echo_only")
                assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                run_phits(Cell(test_surf, material=ins, density=1.0), test_source, control="output_echo_only")
            elif ins.name in ["t-cross", "t-product", "t-time"]:
                inp = make_input([test_cell], test_source, [ins])
                assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                run_phits([test_cell], test_source, [ins]) # TODO: do we want to make this work with output_echo_only?

            # almost everything in misc must have a cell superobject
            elif cls.__module__ == "misc" and cls.name not in ["frag_data", "multiplier"]:
                inp = make_input([Cell([test_surf], test_mat, 0.5, **{cls.name: ins})], test_source, [], control="output_echo_only")
                assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                run_phits([Cell([test_surf], test_mat, 0.5, **{cls.name: ins})], test_source, [], control="output_echo_only")

            elif cls.name == "frag_data":
                inp = make_input([test_cell], test_source, [], control="output_echo_only", cross_sections=[ins])
                assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                run_phits([test_cell], test_source, [], control="output_echo_only", cross_sections=[ins])
            elif cls.name == "multiplier":
                inp = make_input([test_cell], test_source, [], control="output_echo_only", multipliers=[ins])
                assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                run_phits([test_cell], test_source, [], control="output_echo_only", multipliers=[ins])
            elif cls.name == "parameters":
                actual_params = {k: v for k, v in ins.__dict__.items() if k not in ["name", "control", "nuclear_memory_rescale"]}

                inp = make_input([test_cell], test_source, [], control="output_echo_only", **actual_params)
                assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                run_phits([test_cell], test_source, [], control="output_echo_only", **actual_params)
            else:
                raise RuntimeError(f"Class {cls} not handled in test cases.")
        except RuntimeError as e:
            if e.args[0] == "PHITS line limit reached.":
                print(e) # throw away errors that are PHITS's fault
            else:
                raise e

    definition_syntax_correct()


    # The make_input algorithm depends strongly on correctness of the hashing and equality functions defined in base.py
    # @given(builds_right(cls, req, opt), builds_right(cls, req, opt))
    # @settings(deadline=None, suppress_health_check=[HealthCheck.too_slow],
    #           phases=(Phase.explicit, Phase.reuse, Phase.generate, Phase.target,
    #                   Phase.shrink
    #                   ))
    # def eq_is_definitional_equality(ins1, ins2): # this is Noether's first isomorphism theorem
    #     if ins1 == ins2:
    #         assert ins1.definition() == ins2.definition(), "Objects that are __eq__ should have the same definition."

    #     if ins1.definition() == ins2.definition():
    #         assert ins1 == ins2, "Objects that have the same definition should be __eq__."

    # eq_is_definitional_equality()




if __name__ == '__main__':
    import base,  source, parameters, cell, surface, tally, transform, material, misc
    omit = ["Tetrahedral", "MappedMagneticField"] # These require file IO, which shall not be generated.
                                                  # It should however be checked that they do generate Fortran runtime errors,
                                                  # and not parsing errors.
    # omit += ["Gaussian", "Cylindrical", "Rectangular", "GaussianPrism", "Parabolic", "ParabolicPrism", "Spherical", "Beam",
    #          "Conical", "TrianglePrism", "Parameters", "OuterVoid", "Void", "Cell", "MagneticField", "NeutronMagneticField",
    #          "UniformElectromagneticField"]
    # "Tetrahedral", "OuterVoid", "Void", "Cell",
    #         , "MappedMagneticField", "UniformElectromagneticField", "MappedElectromagneticField",
    #         "DeltaRay", "TrackStructure", "SuperMirror", "FragData", "Importance", "WeightWindow", "WWBias", "ForcedCollisions",
    #         "RepeatedCollisions", "Multiplier", "RegionName", "Counter", "Timer", "Plane", "PointPlane", "ParallelPlane",
    #         "Sphere", "Cylinder", "Cone", "SimpleConic", "GeneralConic", "Torus", "Box", "HexagonalPrism",
    #         "EllipticalCylinder", "Spheroid", "Wedge", "TetrahedronBox"]
    # TODO: check Cell correct after fixing Material
    if "Parameters" not in omit:
        test_a(parameters.Parameters)
    for mod in [base, source, parameters, cell, misc, source, surface, tally, transform, material]:
        for name, cls in mod.__dict__.items():
            if name not in omit:
                if isinstance(cls, type) and cls.__module__ == mod.__name__ and issubclass(cls, PhitsObject) and cls != PhitsObject:
                    test_a(cls)
