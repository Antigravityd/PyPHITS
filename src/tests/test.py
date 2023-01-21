parameters = Parameters()
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
