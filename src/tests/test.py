import unittest
from hypothesis import given


def test_run(ob):
    inp = f"""
    [Parameters]
    icntl = 3

    {ob.section_title()}
    {ob.prelude() if hasattr(ob, "prelude") else ''}
    {ob.definition()}
    """



def make_test(cls):
    req = map(lambda tup: tuples(*[i.strat for i in tup[1]]) if isinstance(tup[1], tuple) else tup[1].strat,
              sorted([v for k, v in cls.syntax.items() if v[2] is not None], key=lambda tup: tup[2]))
    opt = map(lambda tup: (tup[0], tuples(*[i.strat for i in tup[1][1]]) if isinstance(tup[1][1], tuple) else tup[1][1].strat),
              [(k, v) for k, v in cls.syntax.items() if v[2] is None])

    @given(builds(cls, *req, **opt))
    def definition_syntax_correct(ins):
        test_source = Cylindrical()
        test_surf = Sphere(1)
        test_cell = Void(test_surf)

        if ins.name == "source":
            run_phits([ins], [test_cell], [], control="input_echo")
        elif ins.name == "cell":
            run_phits([test_source], [ins], [], control="input_echo")
        elif ins.name == "surface":
            run_phits([test_source], [Void(ins)], [], control="input_echo")
        elif ins.name == "material":
            run_phits([test_source], Cell(test_surf, material=ins, density=1.0))





if __name__ == '__main__':
    for
