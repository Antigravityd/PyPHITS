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
    for py_iden, (phits_iden, valspec, idx, *s) in sorted(((k, v) for k, v in cls.syntax.items() if v[2] is not None),
                                                          key=lambda t: t[1][2]):
        # these omitted arguments corrupt the unit-like nature of the tests for no reason
        if py_iden not in ["mask"] and not (hasattr(cls, "subobjects") and py_iden in cls.subobjects):
            if isinstance(valspec, tuple):
                req.append(tuples(*[i.strat for i in valspec]))
            else:
                req.append(valspec.strat)

    opt = dict()
    for py_iden, (phits_iden, valspec, idx, *s) in filter(lambda t: t[1][2] is None, cls.syntax.items()):
        if py_iden not in ["mask"] and not (hasattr(cls, "subobjects") and py_iden in cls.subobjects):
            if isinstance(valspec, tuple):
                opt[py_iden] = one_of(none(), tuples(*[i.strat for i in valspec]))
            else:
                opt[py_iden] = one_of(none(), valspec.strat)



    @composite
    def builds_right(draw, cl, re, op):
        try:
            ob = draw(builds(cl, *re, **op))
        except ValueError as e:
            ob = None

        assume(ob)
        return ob

    @given(builds_right(cls, req, opt))
    @settings(deadline=None, suppress_health_check=[HealthCheck.too_slow, HealthCheck.large_base_example,
                                                    HealthCheck.data_too_large, HealthCheck.filter_too_much],
              phases=(Phase.explicit, Phase.reuse, Phase.generate, Phase.target,
                      # Phase.shrink
                      ))
    def definition_syntax_correct(ins):
        test_source = Cylindrical(["1H"], EnergyDistribution([(1.0, 1.0, 1.0)]))
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
                run_phits([ins, test_cell], test_source, [], control="output_echo_only", automatic_e_bounds=False)

            elif ins.name == "surface":
                inp = make_input([Cell([ins], test_mat, 0.5)], test_source, [], control="output_echo_only")
                assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                run_phits([Cell([ins], test_mat, 0.5)], test_source, [], control="output_echo_only")

            elif ins.name == "material":
                inp = make_input(Cell([test_surf], ins, 0.5), test_source, [], control="output_echo_only")
                assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                run_phits(Cell([test_surf], ins, 0.5), test_source, [], control="output_echo_only")

            elif ins.name in ["t-cross", "t-product", "t-time"]:
                inp = make_input([test_cell], test_source, [ins], control="output_echo_only")
                assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                run_phits([test_cell], test_source, [ins], control="output_echo_only")

            # almost everything in misc must have a cell superobject
            elif (cls.__module__ == "misc" and cls.name not in ["frag_data", "multiplier"]) or cls.__module__ == "transform":
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

            elif cls.name in ["data_max", "mat_time_change", "mat_name_color"]:
                mat = Material([("H", 1)], **{cls.name: ins})
                inp = make_input([Cell([test_surf], mat, 0.5)], test_source, [], control="output_echo_only")
                assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                run_phits([Cell([test_surf], mat, 0.5)], test_source, [], control="output_echo_only")

            elif cls.name == "parameters":
                actual_params = {k: v for k, v in ins.__dict__.items() if k not in ["name", "control", "nuclear_memory_rescale"]}

                inp = make_input([test_cell], test_source, [], control="output_echo_only", **actual_params)
                assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                run_phits([test_cell], test_source, [], control="output_echo_only", **actual_params)

            else:
                raise RuntimeError(f"Class {cls} not handled in test cases.")

        except (RuntimeError, ValueError) as e:
            if e.args[0] == "PHITS line limit reached.":
                print(e) # throw away errors that are PHITS's fault
            elif "Integration problem:" in e.args[0]:
                pass
            else:
                raise e

    definition_syntax_correct()


    # @given(lists(one_of(builds_right(Cell, req, opt), builds_right(Void, req, opt))),
    #        lists(surface_spec.strat),
    #        lists(one_of()))
    # @settings(deadline=None, suppress_health_check=[HealthCheck.too_slow, HealthCheck.large_base_example,
    #                                                 HealthCheck.data_too_large, HealthCheck.filter_too_much],
    #           phases=(Phase.explicit, Phase.reuse, Phase.generate, Phase.target,
    #                   # Phase.shrink
    #                   ))
    # def integrated(cells, surfaces, tallies)



# TODO: fix bug appearing in Cell with material with two 1H slipping past the check
if __name__ == '__main__':
    import base,  source, parameters, cell, surface, tally, transform, material, misc
    omit = ["Parameters"]
    # These require file IO or suffer from poor error reporting/documentation so that I can't implement a Python check for the PHITS error.
    # It should however be checked that they generate only the expected errors.
    omit += ["Tetrahedral", "MappedMagneticField", "MappedElectromagneticField", "FragData"]
    # Takes forever because there's a billion parameters, and they all interdepend.
    # Enable it only to check that things are grammatically correct; perhaps eventually it'll be seamless.
    # omit += ["Parameters", "Gaussian", "Cylindrical", "Rectangular", "GaussianPrism", "Parabolic", "ParabolicPrism",
    #          "Spherical", "Beam", "Conical", "TrianglePrism", "Void", "OuterVoid", "Cell", "MagneticField", "NeutronMagneticField",
    #          "UniformElectromagneticField", "DeltaRay", "TrackStructure", "ElasticOption", "Importance", "WeightWindow", "WWBias",
    #          "ForcedCollisions", "Multiplier", "RegionName", "Counter", "Timer", "Plane", "PointPlane", "ParallelPlane",
    #          "Sphere", "Cylinder", "Cone", "SimpleConic", "GeneralConic", "Box", "EllipticalCylinder", "Spheroid", "Wedge",
    #          "TetrahedronBox", "Transform", "MatNameColor", "Material",  "DumpFluence", "RepeatedCollisions", "DataMax"
    #          ]
    # omit += ["Gaussian", "Cylindrical", "Rectangular", "GaussianPrism", "Parabolic", "ParabolicPrism", "Spherical", "Beam",
    #          "Conical", "TrianglePrism",   "MagneticField", "NeutronMagneticField",
    #          "Void", "OuterVoid", "Cell",
    #          "MagneticField", "NeutronMagneticField", "UniformElectromagneticField",
    #          "DeltaRay", "TrackStructure", "ElasticOption", "Importance", "WeightWindow", "WWBias", "ForcedCollisions",
    #          "Multiplier", "RegionName", "Counter", "Timer", "Plane", "PointPlane", "ParallelPlane", "Sphere", "Cylinder",
    #          "Cone", "SimpleConic", "GeneralConic", "Box", "EllipticalCylinder", "Spheroid", "Wedge", "TetrahedronBox",
    #          "Transform", "DataMax", "MatNameColor", "DumpFluence" # "Material",

    #          "UniformElectromagneticField", "DeltaRay", "TrackStructure", "Importance", "WeightWindow", "WWBias",
    #          "ForcedCollisions", "Multiplier", "RegionName", "Counter", "Timer", "Plane", "PointPlane", "ParallelPlane",
    #          "Sphere", "Cylinder", "Cone", "SimpleConic", "GeneralConic", "Box", "HexagonalPrism", "EllipticalCylinder",
    #          "Spheroid", "Wedge", "TetrahedronBox", "Void", "OuterVoid", # "Cell", # "DumpFluence", "DumpProduct", "DumpTime",
    #          "Transform", "Material"]

    #         "DeltaRay", "TrackStructure", "SuperMirror", "FragData", "Importance", "WeightWindow", "WWBias", "ForcedCollisions",
    #         "RepeatedCollisions", "Multiplier", "RegionName", "Counter", "Timer", "Plane", "PointPlane", "ParallelPlane",
    #         "Sphere", "Cylinder", "Cone", "SimpleConic", "GeneralConic", "Torus", "Box", "HexagonalPrism",
    #         "EllipticalCylinder", "Spheroid", "Wedge", "TetrahedronBox"]
              # ]
    # TODO: check Cell correct after fixing Material
    if "Parameters" not in omit:
        test_a(parameters.Parameters)
    for mod in [base, source, parameters, cell, misc, source, surface, tally, transform, material]:
        for name, cls in mod.__dict__.items():
            if name not in omit:
                if isinstance(cls, type) and cls.__module__ == mod.__name__ and issubclass(cls, PhitsObject) and cls != PhitsObject:
                    test_a(cls)
