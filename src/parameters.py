"""A dumping ground for top-level control parameters."""


from base import *
import sys



class Parameters():
    """A "dictionary with an attitude" representing an entry in the [Parameters] section of an input file.
    Any extra keyword arguments to any constructors are minted into parameter objects.


    >>> print(Parameters(ndedx=2, dbcutoff=3.3).definition())
    ndedx = 2
    dbcutoff = 3.3
    """
    name = "parameters"
    syntax = {"control": ("icntl", FinBij({"normal": 0, "output_cross-section": 1, "output_echo_only": 3, "all_reg_void": 5,
                                           "source_check": 6, "show_geometry": 7, "show_geometry_with_xyz": 8, "show_regions": 9,
                                           "show_regions_with_tally": 10, "show_3d_geometry": 11, "use_dumpall": 12, "sum_tally": 13,
                                           "auto_volume": 14, "ww_bias_tally": 15, "analysis_script": 16, "anatally": 17}), None),
              "max_histories": ("maxcas", PosInt(), None),
              "max_batches": ("maxbch", PosInt(), None),
              "nuclear_memory_rescale": ("xsmemory", PosReal(), None),
              "timeout": ("timeout", NegDisable(), None),
              "stdev_control": ("istdev", FinBij({"history_restart": -2, "batch_restart": -1, "normal": 0, "batch": 1, "history": 2}), None),
              "share_tallies": ("italsh", Choice10(), None),
              "check_consistency": ("ireschk", Choice10(c_style=True), None),
              "xor_prng": ("nrandgen", Choice10(), None),
              # "seed_skip": ("irskeep", Integer(), None),
              "random_seed": ("rseed", Real(), None),
              "seed_from_time": ("itimrand", Choice10(), None),
              # bitrseed?,

              "proton_e_cutoff": ("emin(1)", PosReal(), None),
              "neutron_e_cutoff": ("emin(2)", PosReal(), None),
              "pionp_e_cutoff": ("emin(3)", PosReal(), None),
              "pion0_e_cutoff": ("emin(4)", PosReal(), None),
              "pionm_e_cutoff": ("emin(5)", PosReal(), None),
              "muonp_e_cutoff": ("emin(6)", PosReal(), None),
              "muonm_e_cutoff": ("emin(7)", PosReal(), None),
              "kaonp_e_cutoff": ("emin(8)", PosReal(), None),
              "kaon0_e_cutoff": ("emin(9)", PosReal(), None),
              "kaonm_e_cutoff": ("emin(10)", PosReal(), None),
              "other_e_cutoff": ("emin(11)", PosReal(), None),
              "electron_e_cutoff": ("emin(12)", PosReal(), None),
              "positron_e_cutoff": ("emin(13)", PosReal(), None),
              "photon_e_cutoff": ("emin(14)", PosReal(), None),
              "deuteron_e_cutoff": ("emin(15)", PosReal(), None),
              "triton_e_cutoff": ("emin(16)", PosReal(), None),
              "he3_e_cutoff": ("emin(17)", PosReal(), None),
              "he4_e_cutoff": ("emin(18)", PosReal(), None),
              "nucleon_e_cutoff": ("emin(19)", PosReal(), None),
              "proton_e_max": ("dmax(1)", PosReal(), None),
              "neutron_e_max": ("dmax(2)", PosReal(), None),
              "electron_e_max": ("dmax(12)", PosReal(), None),
              "positron_e_max": ("dmax(13)", PosReal(), None),
              "photon_e_max": ("dmax(14)", PosReal(), None),
              "deuteron_e_max": ("dmax(15)", PosReal(), None),
              "he4_e_max": ("dmax(18)", PosReal(), None),
              # "photonuclear_e_max": ("dpnmax", PosReal(), None),
              # lib(i)??
                  "proton_react_cutoff": ("cmin(1)", PosReal(), None),
              "neutron_react_cutoff": ("cmin(2)", PosReal(), None),
              "pionp_react_cutoff": ("cmin(3)", PosReal(), None),
              "pion0_react_cutoff": ("cmin(4)", PosReal(), None),
              "pionm_react_cutoff": ("cmin(5)", PosReal(), None),
              "muonp_react_cutoff": ("cmin(6)", PosReal(), None),
              "muonm_react_cutoff": ("cmin(7)", PosReal(), None),
              "kaonp_react_cutoff": ("cmin(8)", PosReal(), None),
              "kaon0_react_cutoff": ("cmin(9)", PosReal(), None),
              "kaonm_react_cutoff": ("cmin(10)", PosReal(), None),
              "other_react_cutoff": ("cmin(11)", PosReal(), None),
              "electron_react_cutoff": ("cmin(12)", PosReal(), None),
              "positron_react_cutoff": ("cmin(13)", PosReal(), None),
              "photon_react_cutoff": ("cmin(14)", PosReal(), None),
              "deuteron_react_cutoff": ("cmin(15)", PosReal(), None),
              "triton_react_cutoff": ("cmin(16)", PosReal(), None),
              "he3_react_cutoff": ("cmin(17)", PosReal(), None),
              "he4_react_cutoff": ("cmin(18)", PosReal(), None),
              "nucleon_react_cutoff": ("cmin(19)", PosReal(), None),
              "charged_e_min": ("esmin", PosReal(), None),
              "charged_e_max": ("esmax", PosReal(), None),
              "electron_positron_track_structure_e_min": ("etsmin", PosReal(), None),
              "electron_positron_track_structure_e_max": ("etsmax", PosReal(), None),
              "nucleon_track_structure_e_max": ("tsmax", PosReal(), None),
              "electric_transport_type": ("negs", FinBij({"PHITS": -1, "ignore": 0, "EGS5": 1}), None),
              "automatic_e_bounds": ("nucdata", Choice10(), None),
              "electron_positron_adjust_weight_over_e_max": ("ieleh", Choice10(), None),
              "nucleon_nucleus_model_switch_e": ("ejamnu", PosReal(), None),
              "pion_nucleus_model_switch_e": ("ejampi", PosReal(), None),
              "isobar_max_e": ("eisobar", PosReal(), None),
              "isobar_model": ("isobar", Choice10(), None),
              "bertini_jqmd_switch_e": ("eqmdnu", PosReal(), None),
              "jqmd_e_min": ("eqmdmin", PosReal(), None),
              "jqmd_jamqmd_switch_e": ("ejamqmd", PosReal(), None),
              "incl_control":  ("inclg", FinBij({None: 0, "all": 1, "no_He": 2}), None),
              "icnl_e_min": ("einclmin", PosReal(), None),
              "icnl_e_max": ("einclmax", PosReal(), None),
              "inc_elf_control": ("incelf", Choice10(), None),
              "icnl_elf_e_min": ("eielfmin", PosReal(), None),
              "icnl_elf_e_max": ("eielfmax", PosReal(), None),
              "jqmd_2": ("irqmd", Choice10(), None),
              "scinful_qmd": ("iscinful", Choice10(), None),
              "kerma_mode": ("kerma", Choice10(), None),
              "pseudo_reaction_e": ("epseudo", PosReal(), None),
              "proton_time_cutoff": ("tmax(1)", PosReal(), None),
              "neutron_time_cutoff": ("tmax(2)", PosReal(), None),
              "pionp_time_cutoff": ("tmax(3)", PosReal(), None),
              "pion0_time_cutoff": ("tmax(4)", PosReal(), None),
              "pionm_time_cutoff": ("tmax(5)", PosReal(), None),
              "muonp_time_cutoff": ("tmax(6)", PosReal(), None),
              "muonm_time_cutoff": ("tmax(7)", PosReal(), None),
              "kaonp_time_cutoff": ("tmax(8)", PosReal(), None),
              "kaon0_time_cutoff": ("tmax(9)", PosReal(), None),
              "kaonm_time_cutoff": ("tmax(10)", PosReal(), None),
              "other_time_cutoff": ("tmax(11)", PosReal(), None),
              "electron_time_cutoff": ("tmax(12)", PosReal(), None),
              "positron_time_cutoff": ("tmax(13)", PosReal(), None),
              "photon_time_cutoff": ("tmax(14)", PosReal(), None),
              "deuteron_time_cutoff": ("tmax(15)", PosReal(), None),
              "triton_time_cutoff": ("tmax(16)", PosReal(), None),
              "he3_time_cutoff": ("tmax(17)", PosReal(), None),
              "he4_time_cutoff": ("tmax(18)", PosReal(), None),
              "nucleon_time_cutoff": ("tmax(19)", PosReal(), None),
              "proton_weight_min": ("wc1(1)", PosReal(), None),
              "neutron_weight_min": ("wc1(2)", PosReal(), None),
              "pionp_weight_min": ("wc1(3)", PosReal(), None),
              "pion0_weight_min": ("wc1(4)", PosReal(), None),
              "pionm_weight_min": ("wc1(5)", PosReal(), None),
              "muonp_weight_min": ("wc1(6)", PosReal(), None),
              "muonm_weight_min": ("wc1(7)", PosReal(), None),
              "kaonp_weight_min": ("wc1(8)", PosReal(), None),
              "kaon0_weight_min": ("wc1(9)", PosReal(), None),
              "kaonm_weight_min": ("wc1(10)", PosReal(), None),
              "other_weight_min": ("wc1(11)", PosReal(), None),
              "electron_weight_min": ("wc1(12)", PosReal(), None),
              "positron_weight_min": ("wc1(13)", PosReal(), None),
              "photon_weight_min": ("wc1(14)", PosReal(), None),
              "deuteron_weight_min": ("wc1(15)", PosReal(), None),
              "triton_weight_min": ("wc1(16)", PosReal(), None),
              "he3_weight_min": ("wc1(17)", PosReal(), None),
              "he4_weight_min": ("wc1(18)", PosReal(), None),
              "nucleon_weight_min": ("wc1(19)", PosReal(), None),
              "proton_weight_cutoff": ("wc2(1)", PosReal(), None),
              "neutron_weight_cutoff": ("wc2(2)", PosReal(), None),
              "pionp_weight_cutoff": ("wc2(3)", PosReal(), None),
              "pion0_weight_cutoff": ("wc2(4)", PosReal(), None),
              "pionm_weight_cutoff": ("wc2(5)", PosReal(), None),
              "muonp_weight_cutoff": ("wc2(6)", PosReal(), None),
              "muonm_weight_cutoff": ("wc2(7)", PosReal(), None),
              "kaonp_weight_cutoff": ("wc2(8)", PosReal(), None),
              "kaon0_weight_cutoff": ("wc2(9)", PosReal(), None),
              "kaonm_weight_cutoff": ("wc2(10)", PosReal(), None),
              "other_weight_cutoff": ("wc2(11)", PosReal(), None),
              "electron_weight_cutoff": ("wc2(12)", PosReal(), None),
              "positron_weight_cutoff": ("wc2(13)", PosReal(), None),
              "photon_weight_cutoff": ("wc2(14)", PosReal(), None),
              "deuteron_weight_cutoff": ("wc2(15)", PosReal(), None),
              "triton_weight_cutoff": ("wc2(16)", PosReal(), None),
              "he3_weight_cutoff": ("wc2(17)", PosReal(), None),
              "he4_weight_cutoff": ("wc2(18)", PosReal(), None),
              "nucleon_weight_cutoff": ("wc2(19)", PosReal(), None),
              # "proton_source_weight_min": ("swc1(1)", PosReal(), None),
              # "neutron_source_weight_min": ("swc1(2)", PosReal(), None),
              # "pionp_source_weight_min": ("swc1(3)", PosReal(), None),
              # "pion0_source_weight_min": ("swc1(4)", PosReal(), None),
              # "pionm_source_weight_min": ("swc1(5)", PosReal(), None),
              # "muonp_source_weight_min": ("swc1(6)", PosReal(), None),
              # "muonm_source_weight_min": ("swc1(7)", PosReal(), None),
              # "kaonp_source_weight_min": ("swc1(8)", PosReal(), None),
              # "kaon0_source_weight_min": ("swc1(9)", PosReal(), None),
              # "kaonm_source_weight_min": ("swc1(10)", PosReal(), None),
              # "other_source_weight_min": ("swc1(11)", PosReal(), None),
              # "electron_source_weight_min": ("swc1(12)", PosReal(), None),
              # "positron_source_weight_min": ("swc1(13)", PosReal(), None),
              # "photon_source_weight_min": ("swc1(14)", PosReal(), None),
              # "deuteron_source_weight_min": ("swc1(15)", PosReal(), None),
              # "triton_source_weight_min": ("swc1(16)", PosReal(), None),
              # "he3_source_weight_min": ("swc1(17)", PosReal(), None),
              # "he4_source_weight_min": ("swc1(18)", PosReal(), None),
              # "nucleon_source_weight_min": ("swc1(19)", PosReal(), None),
              "weight_window_max": ("wupn", PosReal(), None),
              "survival_weight": ("wsurvn", PosReal(), None),
              # "max_split": ("mxwpln", PosReal(), None),
              # "window_at": ("mwhere", FinBij({"reaction": -1, "both": 0, "reg_crossing": 1}), None),
              "ww_bias": ("iwwbias", Choice10(), None),
              "std_cutoff": ("istdcut", Choice10(), None),
              # "only_cut_after_batch": ("istdbat", PosInt(), None),
              "stopping_model": ("ndedx", FinBij({"SPAR_nucleus_only+NTMC": 0, "ATIMA+NTMC": 1, "SPAR+NTMC": 2, "ATIMA": 3}), None),
              "atima_db_max": ("mdbatima", PosInt(), None),
              "atima_e_cutoff": ("dbcutoff", PosReal(), None),
              "atima_water_ion_e": ("ih2o", NegDisable(), None),
              "atima_effective_charge": ("ifixchg", Choice10(), None),
              "restricted_delta_LET": ("irlet", Choice10(), None),
              "elastic_scattering": ("ielas", FinBij({None: 0, "neutron": 1, "both": 2}), None),
              "elastic_angle_groups": ("ielms", PosInt(), None),
              "nucleon_model": ("icxnp", FinBij({"JAM": 0, "JENDL": 1}), None),
              "nucleon_nucleus_model": ("icxsni", FinBij({"Perlstein-Niita": 0, "KUROTAMA": 1, "Sato": 2}), None),
              "nucleus_nucleus_model": ("icrhi", FinBij({"Shen": 0, "NASA": 1, "KUROTAMA": 2}), None),
              "deuteron_model": ("icrdm", FinBij({"nuclear": 0, "MWO": 1}), None),
              "pion_induced_model": ("icxspi", FinBij({"geometrical": 0, "Hashimoto": 1}), None),
              "evap_model": ("nevap", FinBij({None: 0, "GEM": 3}), None),
              "gem_version": ("ngem", FinBij({1: 1, 2: 2}), None),
              "fission_model": ("ifission", FinBij({"PHITS": 1, "Iwamoto": 2}), None),
              "gamma_decay_model": ("igamma", FinBij({None: 0, "PHITS": 1, "EBITEM": 2, "EBITEM isomer": 3,
                                                      "PHITS no Doppler": -1, "EBITEM no Doppler": -2, "EBITEM no Doppler": -3}), None),
              "statistical_multifrag_model": ("ismm", Choice10(), None),
              "event_generator": ("e-mode", FinBij({None: 0, "simple": 1, "complicated": 2}), None),
              "event_generator_max_e": ("em-emode", PosReal(), None),
              "neutron_kerma": ("ikerman", Choice10(c_style=True), None),
              "photon_kerma": ("ikermap", Choice10(c_style=True), None),
              "photon_adjoint": ("iadjoint", Choice10(), None),
              "bertini_nucleon_cross_section": ("inmed", FinBij({"free": 0, "cugnon_old": 1, "cugnon_new": 2}), None),
              "bertini_angular_distribution": ("andit",  FinBij({"split": 0, "isotropic": 1, "forward": 2}), None),
              "experimental_neutron_fission": ("iidfs", Choice10(c_style=True), None),
              "discrete_dwba_spectra": ("idwba", Choice10(), None),
              # "absorb_low_e_neg": ("npdik", Choice10(c_style=True), None),
              "neutron_capture_cutoff": ("emcnf", PosReal(), None),
              "fission_delayed_neutrons": ("dnb", OneOf(FinBij({"natural": -1, None: 0}), PosReal()), None),
              "fission_neutron_production": ("nonu", Choice10(), None),
              "s_matrix_interpolation": ("isaba", Choice10(), None),
              "detailed_photon_cutoff": ("emcpf", PosReal(), None),
              "photon_coherent_scatter": ("nocoh", Choice10(c_style=True), None),
              "simple_brem": ("ibad", Choice10(), None),
              "brem_photon_count": ("bnum", PosReal(), None),
              "xray_photon_count": ("xnum", PosReal(), None),
              "substep_brem": ("numb", Choice10(), None),
              # ipegs
              "multiple_scattering": ("imsegs", FinBij({"PHITS-EGS5": 1, "EGS5": 0}), None),
              # iegsout
              "egs_rand_seed": ("iegsrand", PosReal(), None),
              "edge_photons": ("iedgfl", Choice10(), None),
              "edge_auger_electrons": ("iauger", Choice10(), None),
              "coherent_reyleigh_scattering": ("iraylr", Choice10(), None),
              # lpolar
              # iunrst
              "egs_material_size": ("chard", Real(), None),
              "icru90_corrections": ("epstfl", Choice10(), None),
              "densest_gas": ("gasegs", PosReal(), None),
              "compton_incoherent_scattering": ("incohr", Choice10(), None),
              "compton_doppler_broadening": ("iprofr", Choice10(), None),
              "electron_impact_ionization": ("impacr", Choice10(), None),
              "electron_impact_xray_split": ("ieispl", Choice10(), None),
              "electron_impact_xray_count": ("neispl", PosInt(), None),
              "brem_sample_polar_angle": ("ibrdst", Choice10(), None),
              "electron_pair_sample_polar_angle": ("iprdst", Choice10(), None),
              "photoelectron_sample_angle": ("iphter", Choice10(), None),
              "compton_bound_cross_section": ("ibound", Choice10(), None),
              "brem_cross_section_correction": ("iaprim", FinBij({None: 2, "modeled": 1, "empirical": 0}), None),
              "photonuclear_reactions": ("ipnint", FinBij({None: 0, "no_flourescence": 1, "all": 2}), None),
              "photonuclear_probability_weight": ("pnimul", PosReal(), None),
              "photon_induced_muon_production": ("igmuppd", Choice10(), None),
              "muon_capture": ("imucap", FinBij({False: 0, True: 1, "custom_xray": 2}), None),
              "muon_induced_reactions": ("imuint", Choice10(), None),
              "muon_induced_brem": ("imubrm", Choice10(), None),
              "muon_induced_pairs": ("imuppd", Choice10(), None),
              "muon_induced_reaction_min": ("emumin", PosReal(), None),
              "muon_induced_reaction_max": ("emumax", PosReal(), None),
              "neutrino_induced_reaction": ("ntrnore", Choice10(), None),
              "angle_straggling": ("nspred", FinBij({"Lynch": 2, "NTMCC": 1, None: 0}), None),
              "lynch_params": (("ascat1", "ascat2"), (Real(), Real()), None),
              "energy_straggling": ("nedisp", Choice10(), None),
              "gravity": (("gravx", "gravy", "gravz"), (Real(), Real(), Real()), None),
              # usrmgt, usrelst
              "magnetic_field": ("imagnf", Choice10(), None),
              "electromagnetic_field": ("ielctf", Choice10(), None),
              # TODO: decide if users can recover output files from work directory (and, accordingly, if output options are necessary)
              # "binary_cells":   ("icells", FinBij({False: 3, True: }))
              }

    def restrictions(self):
        wmax = self.weight_window_max if self.weight_window_max is not None else 5
        wgt = self.survival_weight if self.survival_weight is not None else 0.6 * wmax
        if not (1 < wgt < wmax):
            raise ValueError(f"One must have 1 < survival_weight < weight_window_max; got survival_weight={wgt}"
                             f"and weight_window_max={wmax}.")

        mhist = self.max_histories if self.max_histories is not None else 10
        mbch = self.max_batches if self.max_batches is not None else 10
        if mhist * mbch > 2_147_483_647:
            raise ValueError(f"History and batch maxima ({self.max_histories} and {self.max_batches}) will result in overflow.")

        for wc1, wc2 in zip(sorted(filter(lambda x: "wc1" in x[1][0], self.syntax.items()),
                                   key=lambda x: re.search("[1-9][0-9]*", x[1][0])[0]),
                            sorted(filter(lambda x: "wc2" in x[1][0], self.syntax.items()),
                                   key=lambda x: re.search("[1-9][0-9]*", x[1][0])[0])):
            v1 = getattr(self, wc1[0])
            v2 = getattr(self, wc2[0])
            if (v1 is not None and v2 is not None and v1 <= v2) or (v1 is None and v2 is not None and v2 >= 0.5):
                raise ValueError(f"One must have {wc1[1][0]} > {wc2[1][0]}; got {wc1[1][0]}={v1} and {wc2[1][0]}={v1}.")


    def __init__(self, **kwargs):
        extra = {k: kwargs[k] for k in set(kwargs) - set(self.syntax)}

        assert not extra, f"Unrecognized parameters {extra} in initialization of Parameters object." \
            "Check that the correct parameters were passed to other objects."
        for k, (_, spec, _) in self.syntax:
            if k in kwargs:
                if isinstance(spec, tuple):
                    map(lambda x: spec.phits(x), kwargs[k])
                else:
                    spec.phits(kwargs[k])

                setattr(self, k, kwargs[k])
            else:
                setattr(self, k, None)

        self.restrictions()

    def __getitem__(self, key):
        return self.__dict__[key]

    def empty(self):
        return True if self.__dict__ == {"name": "parameters"} else False

    def definition(self):
        inp = ""
        for var, val in  self.__dict__.items():
            if val is not None:
                if var in self.syntax:
                    phits_iden, valspec = self.syntax[var][0], self.syntax[var][1]
                    if isinstance(valspec, tuple):
                        for idx, spec in enumerate(valspec):
                            mapped = spec.phits(val[idx])
                            if callable(mapped):
                                raise mapped(var)
                            else:
                                inp += f"{phits_iden[idx]} = {mapped}\n"
                    else:
                        mapped = valspec.phits(val)
                        if callable(mapped):
                            raise mapped(var)
                        else:
                            inp += f"{phits_iden} = {mapped}\n"
        return inp

    def section_title(self) -> str:
        """Return the section title under which a PhitsObject belongs."""
        sec_name = self.name.replace("_", " ").title()
        return f"[{sec_name}]\n"


    @classmethod
    def syntax_desc(self) -> str:
        """Return a readable summary of the initialization syntax of the PhitsObject in question.
        Used to generate documentation, but is useful in interactive sessions to """
        required = sorted([(k, v) for k, v in self.syntax.items() if v[2] is not None], key=lambda tup: tup[1][2])
        opt = [(k, v) for k, v in self.syntax.items() if v[2] is None]
        r = ""
        def capfirst(st):
            return st[0].upper() + st[1:]

        if required:
            r = "Required arguments:\n\n|Position|Python name|PHITS name|Accepted value|\n|----|----|----|----|\n"
            for py_attr, (phits_attr, valspec, position, *s) in required:

                if isinstance(valspec, tuple):
                    j = ", "
                    r += f"|{position}|`{py_attr}`|`{phits_attr}`|A tuple ({j.join(map(lambda x: x.description(), valspec))}).|\n"
                else:

                    r += f"|{position}|`{py_attr}`|`{phits_attr}`|{capfirst(valspec.description())}.|\n"

        if opt:
            r += "\nOptional arguments:\n\n|Python name|PHITS name|Accepted value|\n|----|----|----|\n"
            for py_attr, (phits_attr, valspec, position, *s) in opt:
                if isinstance(valspec, tuple):
                    j = ", "
                    r += f"|`{py_attr}`|`{phits_attr}`|A tuple ({j.join(map(lambda x: x.description(), valspec))}).|\n"
                else:

                    r += f"|`{py_attr}`|`{phits_attr}`|{capfirst(valspec.description())}.|\n"
        return r

    @classmethod
    def syntax_for(self, attr: str) -> str:
        """Return a readable summary of a specific initialization parameter of the PhitsObject in question."""
        r = "PHITS name\tAccepted value\tPosition\n"
        r += f"{self.syntax[attr][0]}\t{self.syntax[attr][1]}\t{self.syntax[attr][2]}"
        return r


__pdoc__ = dict()
__pdoc__["builds"] = False
__pdoc__["slices"] = False
__pdoc__["Parameters"] = Parameters.__doc__ + Parameters.syntax_desc()
