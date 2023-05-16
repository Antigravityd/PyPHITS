import sys
import os
import re

from hypothesis import given, settings, assume, HealthCheck, Phase
from hypothesis.strategies import *

from pyphits import run_phits, make_input, PhitsObject, Cell, Void, Cylindrical, Sphere, Material, EnergyDistribution, \
    FragData, SuperMirror, MatTimeChange, _surface_spec, _source_spec, _cell_spec, _tally_spec
from pyphits.valspec import builds_right, IsA






def test_unitlike():
    def test_a(cls):
        print(cls)
        @given(builds_right(cls, omit=lambda x: x[0] != "mask" and not (hasattr(cls, "subobjects") and x[0] in cls.subobjects)))
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
                    if cls.name == "TetrahedralSource":
                        inp = make_input([test_cell], ins, [], control="output_echo_only")
                        assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                        run_phits([test_cell], ins, [], control="output_echo_only", injected_files={ins.tetreg.tet_file: "Bruh"})
                    else:
                        inp = make_input([test_cell], ins, [], control="output_echo_only")
                        assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                        run_phits([test_cell], ins, [], control="output_echo_only")

                elif ins.name == "cell":
                    if cls.name == "Tetrahedral":
                        inp = make_input([ins, test_cell], test_source, [], control="output_echo_only")
                        assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                        run_phits([ins, test_cell], test_source, [], control="output_echo_only", injected_files={ins.tet_file: "Bruh"})
                    else:
                        inp = make_input([ins, test_cell], test_source, [], control="output_echo_only")
                        assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                        run_phits([ins, test_cell], test_source, [], control="output_echo_only")

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


                elif hasattr(cls, "superobjects") and "cell" in cls.superobjects or cls.name == "transform":
                    if cls.__name__ == "MappedMagneticField":
                        inp = make_input([Cell([test_surf], test_mat, 0.5, **{cls.name: ins})], test_source, [],
                                         control="output_echo_only")
                        assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                        run_phits([Cell([test_surf], test_mat, 0.5, **{cls.name: ins})], test_source, [],
                                         control="output_echo_only", injected_files={ins.m_file: "Bruh"})

                    elif cls.__name__ == "MappedElectromagneticField":
                        inp = make_input([Cell([test_surf], test_mat, 0.5, **{cls.name: ins})], test_source, [],
                                         control="output_echo_only")
                        assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                        run_phits([Cell([test_surf], test_mat, 0.5, **{cls.name: ins})], test_source, [],
                                         control="output_echo_only", injected_files={ins.e_file: "Bruh", ins.m_file: "Bruh"})
                    else:
                        inp = make_input([Cell([test_surf], test_mat, 0.5, **{cls.name:ins})], test_source, [], control="output_echo_only")
                        assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                        run_phits([Cell([test_surf], test_mat, 0.5, **{cls.name: ins})], test_source, [], control="output_echo_only")


                elif cls.name in ["frag_data", "multiplier", "super_mirror", "mat_time_change"]:
                    if cls.name == "frag_data":
                        inp = make_input([test_cell], test_source, [], control="output_echo_only", **{cls.name + "s": ins})
                        assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                        run_phits([test_cell], test_source, [], control="output_echo_only", **{cls.name + "s": ins},
                                         injected_files={ins.file: ""})
                    else:
                        inp = make_input([test_cell], test_source, [], control="output_echo_only", **{cls.name + "s": ins})
                        assert ins.section_title() in inp, f"Section didn't make it into input:\n{inp}"
                        run_phits([test_cell], test_source, [], control="output_echo_only", **{cls.name + "s": ins})


                elif cls.name in ["data_max", "mat_time_change", "mat_name_color",]:
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
                elif "Integration problem:" in e.args[0]: # or that result from global_restrictions successfully finding an error
                    pass
                else:
                    raise e

        definition_syntax_correct()

    import pyphits

    # Parameters is too slow/too much integration dependence, and the distributions are tested integrated with Cells
    omit = ["Parameters", "TimeDistribution", "AngleDistribution", "EnergyDistribution",
            "MappedMagneticField", "MappedElectromagneticField", "FragData", "Tetrahedral", "TetrahedralSource"] # these need file IO
    # Progressively uncomment these during debugging to avoid waiting on redundant tests
    omit += ["Transform", "Plane", "PointPlane", "ParallelPlane", "Sphere", "Cylinder", "Cone", "SimpleConic", "GeneralConic",
             "Torus", "Box", "EllipticalCylinder", "Spheroid", "Wedge", "TetrahedronBox", "MagneticField", "NeutronMagneticField",
             "ElectromagneticField", "DeltaRay", "TrackStructure", "ElasticOption", "Importance", "WeightWindow", "WWBias",
             "ForcedCollisions", "RepeatedCollisions", "RegionName", "Counter", "Timer", "DataMax", "MatNameColor",
             "Material",  "Void", "OuterVoid", "Cell", # "SuperMirror", "MatTimeChange", "Cylindrical", "Rectangular", "Gaussian",
             # "GaussianSlices", "Parabolic",
             # "ParabolicSlices", "Spherical", "Beam", "Conical", "TriangularPrism", "SurfaceSource", # "DumpFluence", "DumpProduction",
             # "DumpTime"
             ]
    if "Parameters" not in omit:
        test_a(parameters.Parameters)
    for name, cls in list(pyphits.__dict__.items()):
        if name not in omit:
            if isinstance(cls, type) and cls.__module__ == pyphits.__name__ and issubclass(cls, PhitsObject) and cls != PhitsObject:
                test_a(cls)

if __name__ == "__main__":
    test_unitlike()
