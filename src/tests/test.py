parameters = Parameters(icntl=7)
sources = [Cylindrical("neutron", 2.2),
           Rectangular("neutron", 12),
           Gaussian("electron", 11),
           Parabolic("proton", 13),
           Spherical("208Pb", 1.5),
           # Beam,
           Conical("Am-241", 125),
           Prisim("U-235", 100)
           # Grid,
           # Tetrahedral,
           # SurfaceSource,
           ]
materials = []
surfaces = []
cells = []
transforms = []
#temperatures = []
mat_time_changes = []
magnetic_fields = []
neutron_magnetic_fields = []
uniform_electromagnetic_field = []
mapped_electromagnetic_field = []
delta_ray = []
track_structure = []
importance = []
wi

import sys
import os
sys.path.append(os.getcwd() + '/..')


from base import *

examples = {"parameters":
            [Parameters(r"""icntl = 7
            emin(2) = 1.3e3
            """)],

            "time_distribution_rectangle":
            [Rectangle(r"""t-type = 1
            tw = 2.3
            tn = 2
            td = 0.2
            t0 = 12
            """),
             Rectangle(r"""t-type = 1
             tn = 2
             td = 0.3
             t0 = -2.3
             """)],

            "time_distribution_gaussian":
            [Gaussian(r"""t-type = 2
            tw = 2.3
            tn = 2
            td = 0.2
            tc = 12
            """),
             Gaussian(r"""t-type = 2
             tw = 2.3
             tn = 2
             td = 0.2
             """)],

            "time_distribution_bins":
            [TimeDistributionBins(r"""t-type = 4
            ntt = 3
            1 2 3
            4 5 6
            o-type = 1
            1 2 3 4 5 6
            """),
             TimeDistributionBins(r"""t-type = 3
             ntt = 3
             1 2 3
             4 5 6
             """)],

            "time_distribution_function":
            [TimeDistributionFunction(r"""t-type = 6
            tg1 = -0.30
            tg2 = 0.30
            h(x) = x ** 2
            ll = 30
            o-type = 1
            1 2 3 4 5 6
            """),
             TimeDistributionFunction(r"""t-type = 6
             tg1 = -0.30
             tg2 = 0.30
             h(x) = x ** 2
             ll = 30
             o-type = 1
             """),
             TimeDistributionFunction(r"""t-type = 5
             tg1 = -0.30
             tg2 = 0.30
             h(x) = x ** 2
             ll = 30
             """)],

            "angle_distribution_function":
            [Angle]}

# how to test???
def grammars_are_valid():
    for typ, lis in examples:
        for po in lis:
            Lark(po.grammar())

def definition_inverts_parse(po):
    pass
