"""This is the documentation for PyPHITS, a third-party Python interface to JAEA's \
[Particle Heavy-Ion Transport Code System (PHITS)](https://phits.jaea.go.jp).


See the [README](https://github.com/Antigravityd/PyPHITS/README.org) for an introduction.
"""

# TODOs:
# Counter/multiplier in tally
# Radioisotope and CosmicRay energy distributions
# Dir = iso for spherical shell sources (necessary for CosmicRay to be used properly)
# readable_remapping
# nuclide-specific library settings

from pyphits.valspec import *

import sys
import numpy as np
import collections as col
from numpy.linalg import det
from datetime import datetime
import subprocess as sp
import itertools as it
from functools import reduce
import pandas as pd
import tempfile as tf
import re
import os
from copy import deepcopy
from fortranformat import FortranRecordReader
from shutil import copy

# Configuration options
readable_remapping = True
"""If `True`, use a remapped initialization syntax that's more informative.
   E.g., `Parameters(control="output_echo_only")` instead of `Parameters(icntl=3)`.

   If `False`, use a syntax as close as possible to that specified in the PHITS manual (i.e., the latter form in the above example).
   Note that some identifiers in PHITS do not conform to Python's identifier syntax (e.g. 2d-type, emin(14));
   these identifiers are sanitized as follows:

       - dashes -> underscores
       - beginning with a number -> that clause moved to the end
       - parentheses -> omitted"""


## PARAMETERS
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
        for k, (_, spec, _) in self.syntax.items():
            if k in kwargs and kwargs[k] is not None:
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

    @classmethod
    def required_args(self):
        """The required arguments from `self.syntax`, in order."""
        return list(sorted(((k, v) for k, v in self.syntax.items() if v[2] is not None), key=lambda t: t[1][2]))

    @classmethod
    def optional_args(self):
        """The optional arguments from `self.syntax`."""
        return list(filter(lambda t: t[1][2] is None, self.syntax.items()))


def _tuplify(xs: list) -> tuple:
    """Convert lists to tuples for hashability purposes."""
    return tuple(map(lambda x: _tuplify(x) if isinstance(x, list) else x, xs))

def _continue_lines(inp: str) -> str:
    """If a line is too long for PHITS to handle, attempt to line-break at whitespace.
    Works only for [Surface] and [Cell]."""
    r = ""
    for line in inp.split("\n"):
        if len(line) > 175:
            words = line.split(" ")
            length = 0

            remain = list(it.takewhile(lambda x: len(x[1]) < 175, enumerate(it.accumulate(words, lambda x, y: x + " " + y))))
            contin = " ".join(words[remain[-1][0]+1:])
            remain = remain[-1][1]
            if contin == "" or contin.isspace():
                r += remain + "\n"
            else:
                r += remain + "\n     " + _continue_lines(contin)
        else:
            r += line + "\n"


    return r



## BASE
class PhitsObject:
    """The base factory class distinguishing objects that are intended to end up in some section of a `.inp` file,
    and defining equality and hashability of such objects in sensible ways.

    PhitsObject values always correspond to some section of an input file.
    """

    name = None
    """A string corresponding to the PHITS .inp section an object appears in. See `PhitsObject.names`"""

    syntax = dict()
    """A dictionary with entries of the form `"python_identifier": ("PHITS_identifier", acceptable_ValSpec, arg_position, Optional(none_val))`.
    The key is the attribute on the Python PhitsObject instance, and the keyword argument necessary to set it.
    The second is what type the object must be, as a ValSpec.
    The third is what index in `*args` the argument must have; if it is set to `None`, the argument is optional.
    The last is a value to put in the .inp if the attribute is `None` at compile-time; ordinarily, it's just nothing.
    The first two arguments can be tuples, in which case the passed value must be a tuple of the specified type;
    the single Python assignment corresponds to the entrywise assignments of this tuple to the PHITS identifiers in the .inp."""

    shape = tuple()
    """A tuple that details how the object is to be represented in the .inp file.
    The syntax is inspired by Emacs Lisp skeletons, with some optimizations to tailor it for this use-case.
    In general, the strings in the tuple are inserted verbatim, with a newline in between them.
    If an entry is the name of an attribute of the instance, say `"attr"` whose value is not `None`,
    then the string `PHITS_identifier = acceptable_ValSpec.phits(self.attr)`, based on a lookup of `attr` in `PhitsObject.syntax`, is inserted.
    If an entry is `self`, then `self.index` is inserted.
    To avoid the above two behaviors, and insert a string verbatim, prepend a quote `'`.
    If an entry is another tuple, that tuple gets evaluated as above, but the `identifier = ` is not inserted,
    and a mere space separates the entries.
    Appending a `\\` to any string disables the insertion of the spacing that would otherwise follow;
    similarly, having `\\` as the last string of a tuple entry disables the newline that would otherwise follow.

    If this attribute is callable, it gets called on the instance of the PhitsObject, and the result is processed according to the rules above.
    """

    index = None
    """The compile-time assigned number of an object. Should not be set directly."""

    # no_hash = {"index", "value_map", "ident_map", "nones", "shape", "subobjects", "required", "positional", "optional",
    #            "group_by", "separator", "prelude", "max_groups", "group_size", "parser", "validator"}
    # """Attributes that don't affect the identity of a PhitsObject."""

    names = {"parameters", "source", "material", "surface", "cell", "transform", "temperature","mat_time_change","magnetic_field",
             "electromagnetic_field", "frag_data", "data_max",
             "delta_ray", "track_structure", "super_mirror", "elastic_option", "importance", "weight_window", "ww_bias",
             "forced_collisions", "repeated_collisions", "volume", "multiplier", "mat_name_color", "reg_name", "counter", "timer",
             "t-track", "t-cross", "t-point", "t-adjoint", "t-deposit", "t-deposit2", "t-heat", "t-yield", "t-product", "t-dpa",
             "t-let", "t-sed", "t-time", "t-interact", "t-dchain", "t-wwg", "t-wwbg", "t-volume", "t-gshow", "t-rshow","t-3dshow"}
    """The permissible `PhitsObject.name`s."""

    group_by = None
    """A key function by which to group PhitsObjects together.
    These groups are necessary, for instance, in the `Importance` section, where one can set different importances for different particles
    in the same cell."""

    separator = None
    """A function returning a string to be placed between the definitions of the groups; see `PhitsObject.group_by`."""

    max_groups = None
    """The maximum number of groups of objects of a given type; see `PhitsObject.group_by`."""

    prelude = tuple()
    """A skeleton just like `PhitsObject.shape`, but inserted before all definitions of the `PhitsObject` subclass in question."""
    subobjects = []
    """A list of all `PhitsObject.name`s that can appear as attributes of the object in question."""
    superobjects = []
    """A list of all `PhitsObject.name`s that this object ought to be defined from."""
    restrictions = lambda self: tuple()
    """A function that's called right at the end of initialization. Should be used to raise errors for any combination of values
    that is incorect but not easily expressible via `ValSpec`s."""
    def __init__(self, *args,  **kwargs):
        """Arguments are interpreted according to `PhitsObject.syntax`, and then any leftovers in `kwargs` are used to create a `Parameters`
        object, which is then assigned to a `parameters` attribute."""
        assert self.name in self.names, f"Unrecognized PHITS type {self.name} in PhitsObject initialization."

        # Handle required args
        required = list(map(lambda tup: tup[0], self.required_args()))

        assert len(args) == len(required), f"Wrong number of positional arguments specified in the definition of {self} object."
        for idx, arg in enumerate(args):
            # Validate first
            details = self.syntax[required[idx]]
            valspec = details[1]
            if isinstance(valspec, tuple):
                for id2, spec in enumerate(valspec):
                    mapped = spec.phits(arg[id2])
                    if callable(mapped):
                        raise mapped(required[idx])
            else:
                mapped = valspec.phits(arg)
                if callable(mapped):
                    raise mapped(required[idx])

            # Then set the attribute
            setattr(self, required[idx], arg if not isinstance(arg, list) else _tuplify(arg))

        # Handle optional args
        for arg in self.syntax:
            valspec, position = self.syntax[arg][1], self.syntax[arg][2]
            if position == None:
                if arg in kwargs:
                    # Validate first
                    if kwargs[arg] is not None:
                        if isinstance(valspec, tuple):
                            for idx, spec in enumerate(valspec):
                                mapped = spec.phits(kwargs[arg][idx])
                                if callable(mapped):
                                    raise mapped(arg)
                        else:
                            mapped = valspec.phits(kwargs[arg])
                            if callable(mapped):
                                raise mapped(arg)

                        # Then set the attribute
                        setattr(self, arg, kwargs[arg] if not isinstance(kwargs[arg], list) else _tuplify(kwargs[arg]))
                    else:
                        setattr(self, arg, None)
                else:
                    # Handle unpassed
                    setattr(self, arg, None)


        for attr in self.subobjects:
            if hasattr(self, attr):
                child = getattr(self, attr)
                if child is not None:
                    setattr(child, self.name, self)

        # for attr in self.superobjects:
        #     setattr(self, attr, None) # subobjects are instantiated before superobjects

        remaining = {k: v for k, v in kwargs.items() if k not in self.syntax and k not in self.subobjects}
        if remaining:
            self.parameters = Parameters(**remaining)

        # check restrictions satisfied
        self.restrictions()

    def _add_definition(self, how: tuple, to: str, assignments: bool = True) -> str:
        """Recursively performs skeleton insertions according to `how` at the end of `to`."""
        if callable(how):
            how = how()

        for attr in how:

            endstr = "\n" if assignments else " "
            spacing = " "

            if isinstance(attr, str) and len(attr) > 0 and attr[0] == "'":
                to += attr[1:]
                to += endstr
                continue

            if isinstance(attr, str) and len(attr) > 0 and attr[-1] == "\\":
                endstr = " "
                spacing = ""
                attr = attr[:-1]

            if isinstance(attr, tuple):
                to += self._add_definition(attr, "", assignments=False)
                if attr[-1] == "\\":
                    to += " "
                else:
                    to += "\n"

            elif attr in self.syntax:
                val = getattr(self, attr)
                phits_iden = self.syntax[attr][0]
                valspec = self.syntax[attr][1]
                noneval = ""
                if len(self.syntax[attr]) > 3:
                    noneval = self.syntax[attr][3]

                if val is not None:
                    if isinstance(phits_iden, tuple):
                        for i, (phits, spec) in enumerate(zip(phits_iden, valspec)):
                            assign = f"{phits}{spacing}={spacing}" if assignments else ""
                            v = spec.phits(val[i])
                            if callable(v):
                                raise v(attr)
                            else:
                                to += f"{assign}{v}{endstr}"
                    else:
                        assign = f"{phits_iden}{spacing}={spacing}" if assignments else ""
                        v = valspec.phits(val)
                        if callable(v):
                            raise v(attr)
                        else:
                            to += f"{assign}{v}{endstr}"

                elif noneval != "" and noneval is not None: # I think we don't use nones for anything except the simplest case
                    to += str(noneval)
                    to += endstr


            else:
                if attr == "self":
                    to += str(self.index)
                    to += endstr

                elif attr in self.superobjects: # can't use syntax because impredicativity in module imports
                    to += str(getattr(self, attr).index)
                    to += endstr
                else:
                    to += attr
                    to += endstr
        return to

    def prelude_str(self) -> str:
        """Return a string to appear before the collection of all definitions of subclass instances in an `.inp` file."""
        inp = self._add_definition(self.prelude, "")

        return inp

    def definition(self) -> str:
        """Return the string representing the particular PhitsObject in an `.inp` file."""
        inp = self._add_definition(self.shape, "")
        if self.name in ["surface", "cell"]:
            return _continue_lines(inp)
        else:
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
    def syntax_for(self, attr: str, phits=False) -> str:
        """Return a readable summary of a specific initialization parameter of the PhitsObject in question."""
        if phits:
            rev = self._syntax_reversed()
            r = "Python name\tAccepted value\tPosition\n"
            r += f"{rev[attr][0]}\t{rev[attr][1]}\t{rev[attr][2]}"
            return r
        else:
            r = "PHITS name\tAccepted value\tPosition\n"
            r += f"{self.syntax[attr][0]}\t{self.syntax[attr][1]}\t{self.syntax[attr][2]}"
            return r

    def _sanitize(self, iden: str) -> str:
        """Make a PHITS identifier Python-acceptable.
        See `readable_syntax`.
        """
        r = iden.replace("-", "_")
        r = r.replace("(", "")
        r = r.replace(")", "")
        if r[0].isdigit():
            if m := re.match(".*?_", r):
                r = r[m.end():] + "_" + r[:m.end() - 1]

        return r

    def _syntax_reversed(self):
        """Produce a dictionary analogous to `PhitsObject.syntax` but using PHITS's names."""
        return {v[0]: (k, *v[1:]) for k, v in self.syntax.items()}

    def __eq__(self, other: "PhitsObject") -> bool:
        """PhitsObjects should be equal if their definitions would be the same."""
        if type(self) != type(other):
            return False
        elif hasattr(self, "__dict__") and hasattr(other, "__dict__"):
            d1 = {k: v for k, v in self.__dict__.items() if k in self.syntax}
            d2 = {k: v for k, v in other.__dict__.items() if k in self.syntax}
            return d1 == d2

    def __hash__(self) -> int:
        """PhitsObjects have the hash of their identity-defining attributes."""
        return hash(tuple(v for k, v in sorted(self.__dict__.items())
                          if k in self.syntax and (self not in v.__dict__.values() if hasattr(v, "__dict__") else True)))

    def _repr_pretty_(self, p, cycle) -> str:
        """Hypothesis uses this when printing failing cases."""
        try:
            p.text(str(self.__dict__) + "\n\n" + self.definition())
        except:
            p.text(str(self.__dict__))

    @classmethod
    def required_args(self):
        """The required arguments from `self.syntax`, in order."""
        return list(sorted(((k, v) for k, v in self.syntax.items() if v[2] is not None), key=lambda t: t[1][2]))

    @classmethod
    def optional_args(self):
        """The optional arguments from `self.syntax`."""
        return list(filter(lambda t: t[1][2] is None, self.syntax.items()))


## DISTRIBUTION
class TimeDistribution(PhitsObject):
    """An arbitrary distribution of source weights over time."""
    name = "source"
    syntax = {"bins": (None, List(Tuple(PosReal(), PosReal())), 0),
              "last_bin": (None, PosReal(), 1),
              "particle_production": (None, List(PosReal()), None)}

    shape = lambda slf: ("t-type = 4" if slf.particle_production else "t-type = 3",
                         f"ntt = {len(slf.bins)}",
                         "\n".join(zastr(j[0]) + " " + str(j[1]) for j in slf.bins),
                         ("last_bin",),
                         "o-type = 1\n" + " ".join((str(i) for i in slf.particle_production)) if slf.particle_production else "")


class AngleDistribution(PhitsObject):
    """An arbitrary distribution of source weights over angle."""
    name = "source"
    syntax = {"bins": (None, List(Tuple(PosReal(), PosReal())), 0),
              "last_bin": (None, PosReal(), 1),
              "unit": ("a-type", FinBij({"cos": 1, "degree": 11}), None),
              "particle_production": (None, List(PosReal()), None)}

    shape = lambda self: ("a-type = 14" if self.unit == "degree" else "a-type = 4",
                          f"na = {len(self.bins)}",
                          "\n".join(" ".join(str(i) for i in j) for j in self.bins) + f"\n{self.last_bin}",
                          "q-type = 1\n" + " ".join((str(i) for i in self.particle_production)) if self.particle_production \
                          else "")

    def restrictions(self):
        if self.particle_production is not None and len(self.bins) != len(self.particle_production):
            raise ValueError("For EnergyDistribution: len(bins) must equal len(particle_production).")



class EnergyDistribution(PhitsObject):
    """An arbitrary distribution of source weights over energy."""
    name = "source"
    syntax = {"bins": (None, List(Tuple(PosReal(), PosReal(), PosReal())), 0),
              "adjust": (None, FinBij({"particles": "particles", "weights": "weights"}), None),
              "units": (None, FinBij({"MeV": "MeV", "Angstrom": "Angstrom"}), None),
              "normalize": (None, FinBij({"1/Lethargy": -1, "1/MeV": 1}), None), # TODO: check 1/MeV
              "particle_production": (None, List(PosReal()), None)
              }

    shape = lambda slf: (("e-type = 22" if slf.units == "MeV" else "e-type = 32") if slf.adjust == "particles" else \
                          ("e-type = 23" if slf.units == "MeV" else "e-type = 33"),
                          f"ne = -{len(slf.bins)}" if slf.normalize == "1/Lethargy" else f"ne = {len(slf.bins)}",
                          "\n".join(" ".join(str(j) for j in i) for i in slf.bins),
                          "p-type = 1\n" + " ".join((str(i) for i in slf.particle_production)) if slf.particle_production else "")
    def restrictions(self):
        if self.particle_production is not None and len(self.bins) != len(self.particle_production):
            raise ValueError("For EnergyDistribution: len(bins) must equal len(particle_production).")

# TODO: radioisotope, cosmicray

## TRANSFORM


class Transform(PhitsObject): #
    """An \\(\\mathbb{R}^3\\) isometry represented by a translation vector and a rotation matrix.
    """
    name = "transform"
    syntax = {"translation": (None, Tuple(Real(), Real(), Real()), 0),
              "rotation": (None, Tuple(Real(), Real(), Real()), 1),
              "rotate_first": (None, FinBij({True: 2, False: -2}), None, -2),
              "units": (None, FinBij({"degrees": "degrees", "radians": "radians"}), None)}
    shape = lambda self: ((f"*TR{self.index}" if self.units == "degrees" else f"TR{self.index}",
                           " ".join(str(i) for i in self.translation),
                           " ".join(str(i) for i in self.rotation),
                           "0 0 0 0 0 0",
                           "rotate_first"),)

## SURFACE

surface_common = {"reflective": (None, Choice10(), None),
          "white": (None, Choice10(), None),
          "transform": (None, IsA(Transform, index=True), None),
          "inside": (None, Choice10(), None)}


class Plane(PhitsObject):
    """A plane of the form Ax + By + Cz - D = 0."""
    name = "surface"
    syntax = surface_common | {"A": (None, Real(), 0),
                       "B": (None, Real(), 1),
                       "C": (None, Real(), 2),
                       "D": (None, Real(), 3)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "P", "A", "B", "C", "D"),)

    def restrictions(self):
        if self.A == 0 and self.B == 0 and self.C == 0:
            raise ValueError("For Plane: at least one of A, B, or C must be nonzero.")


# TODO: consider obliterating the next 2
class PointPlane(PhitsObject):
    """A plane specified by three points."""
    name = "surface"
    syntax = surface_common | {"p1": (None, Tuple(Real(), Real(), Real()), 0),
                       "p2": (None, Tuple(Real(), Real(), Real()), 1),
                       "p3": (None, Tuple(Real(), Real(), Real()), 2)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "P",
                           f"{self.p1[0]}", f"{self.p1[1]}", f"{self.p1[2]}",
                           f"{self.p2[0]}", f"{self.p2[1]}", f"{self.p2[2]}",
                           f"{self.p3[0]}", f"{self.p3[1]}", f"{self.p3[2]}"),)

    def restrictions(self):
        if self.p1[0] * (self.p2[1] - self.p3[1]) + self.p2[0] * (self.p3[1] - self.p1[1]) \
           + self.p3[0] * (self.p1[1] - self.p2[1]) == 0: # i.e. points are colinear
            raise ValueError("For PointPlane: p1, p2, and p3 must not line on a line;"
                             f" got p1={self.p1}, p2={self.p2}, and p3={self.p3}.")


class ParallelPlane(PhitsObject):
    """A plane of the form x_i = D."""
    name = "surface"
    syntax = surface_common | {"parallel": (None, FinBij({"x": "X", "y": "Y", "z":"Z"}), 0),
                       "D": (None, Real(), 1)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", f"P{self.parallel}", "D"),)


class Sphere(PhitsObject):
    "A sphere of radius R centered on (x0, y0, z0)."
    name = "surface"
    syntax = surface_common | {"radius": (None, PosReal(), 0),
                       "center": (None, Tuple(Real(), Real(), Real()), 1)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform",
                           "SPH", f"{self.center[0]}", f"{self.center[1]}", f"{self.center[2]}", "radius"),)


class Cylinder(PhitsObject):
    """A right-circular cylinder with center of the bottom face (x_0, y_0, z_0), height vector from the bottom to top face (H_x, H_y, H_z),
    and radius R."""
    name = "surface"
    syntax = surface_common | {"center": (None, Tuple(Real(), Real(), Real()), 0),
                       "height": (None, Tuple(Real(), Real(), Real()), 1),
                       "radius": (None, PosReal(), 2)}
    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform",
                           "RCC", " ".join(str(i) for i in self.center), " ".join(str(i) for i in self.height), "radius"),)

    def restrictions(self):
        if self.height == (0, 0, 0):
            raise ValueError("Cylinder must have a nonzero height vector.")

class Cone(PhitsObject):
    """A truncated right-angle cone with bottom-face center (x_0, y_0, z_0), height vector (H_x, H_y, H_z), and bottom and top radii
    R_1 and R_2 respectively."""
    name = "surface"
    syntax = surface_common | {"center": (None, Tuple(Real(), Real(), Real()), 0),
                               "height": (None, Tuple(Real(), Real(), Real()), 1),
                               "radii": (("top", "bottom"), Interval(0), 2)}
    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform",
                           "TRC", "center", "height", "radii"),)

    def restrictions(self):
        if self.height == (0, 0, 0):
            raise ValueError("Cone must have a nonzero height vector.")




class SimpleConic(PhitsObject): # ellipsoid, hyperboloid, or paraboloid parallel to an axis of the form
                   # A(x-x0)^2+B(y-y0)^2+C(z-z0)^2+2D(x-x0)+2E(y-y0)+2F(z-z0)+G = 0
    name = "surface"
    syntax = surface_common | {"quadratic": ((None, None, None), (Real(), Real(), Real()), 0),
                       "linear": ((None, None, None), (Real(), Real(), Real()), 1),
                       "constant": (None, Real(), 2),
                       "center": ((None, None, None), (Real(), Real(), Real()), 3)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "SQ", "quadratic", "linear", "constant", "center"),)



class GeneralConic(PhitsObject): # ellipsoid, hyperboloid, or paraboloid of the form
                    # A(x-x0)^2+B(y-y0)^2+C(z-z0)^2+Dxy+Eyz+Fzx+Gx+Hy+Jz+K = 0
    name = "surface"
    syntax = surface_common | {"quadratic": ((None, None, None), (Real(), Real(), Real()), 0),
                       "mixed": ((None, None, None), (Real(), Real(), Real()), 1),
                       "linear": ((None, None, None), (Real(), Real(), Real()), 2),
                       "constant": (None, Real(), 2)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "GQ", "quadratic", "mixed", "linear", "constant"),)

# TODO: I don't know what "skewed" means for transfomations on tori, so disabling for now.
class Torus(PhitsObject): # torus parallel to an axis of the form
             # (axisvar - axis0)^2/B^2 + (quadrature(<non-axis displacements>) - A)^2 - 1 = 0
    name = "surface"
    syntax = surface_common | {"axis": (None, FinBij({"x": "X", "y": "Y", "z":"Z"}), 0),
                       "center": ((None, None, None), (Real(), Real(), Real()), 1),
                       "scales": ((None, None, None), (PosReal(), PosReal(), PosReal()), 2)}
    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           f"T{self.axis}", "center", "scales"),)



class Box(PhitsObject): # box formed by three vectors with tails at a given base point, or cross product of 3 intervals,
           # stored in the form x0 y0 z0 Ax Ay Az Bx By Bz Cx Cy Cz
    name = "surface"
    syntax = surface_common | {"base": ((None, None, None), (Real(), Real(), Real()), 0),
                       "walls": (None, OrthogonalMatrix(), 1)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "BOX", "base", " ".join(" ".join(str(i) for i in j) for j in self.walls)),)


# Too complicated to generate
# class HexagonalPrism(PhitsObject):
#     name = "surface"
#     syntax = surface_common | {"base": ((None, None, None), (Real(), Real(), Real()), 0),
#                        "height": ((None, None, None), (Real(), Real(), Real()), 1),
#                        "s1": ((None, None, None), (Real(), Real(), Real()), 2),
#                        "s2": ((None, None, None), (Real(), Real(), Real()), 3),
#                        "s3": ((None, None, None), (Real(), Real(), Real()), 4)}

#     shape = lambda self: ((f"*{self.index}" if self.reflective else
#                            (f"+{self.index}" if self.white else f"{self.index}"),
#                            "transform", "HEX", "base", "height", "s1", "s2", "s3"),)

#     def restrictions(self):
#         if self.height == (0, 0, 0) or self.s1 == (0, 0, 0) or self.s2 == (0, 0, 0) or self.s3 == (0, 0, 0):
#             raise ValueError("HexagonalPrism must have a nonzero height vector.")


class EllipticalCylinder(PhitsObject):
    name = "surface"
    syntax = surface_common | {"center": ((None, None, None), (Real(), Real(), Real()), 0),
                       "axes": (None, OrthogonalMatrix(), 1)}
    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "REC", "center", " ".join(" ".join(str(i) for i in j) for j in self.axes)),)



class Spheroid(PhitsObject):
    name = "surface"
    syntax = surface_common | {"focus1": ((None, None, None), (Real(), Real(), Real()), 0),
                       "focus2": ((None, None, None), (Real(), Real(), Real()), 1),
                       "major_axis": (None, Real(), 2)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "ELL", "focus1", "focus2", "major_axis"),)
    def restrictions(self):
        if self.focus1 == self.focus2:
            raise ValueError(f"Spheroid must have distinct foci; got focus1={self.focus1} and focus2={self.focus2}.")

        if self.major_axis == 0:
            raise ValueError("Spheroid must have a nonzero major axis length.")

        if self.major_axis - np.linalg.norm(np.array(self.focus1) - np.array(self.focus2)) <= 0:
            raise ValueError("Spheroid must have nonzero major axis length larger than the distance betwen its foci;"
                             f" got major_axis={self.major_axis}, focus1={self.focus1}, and focus2={self.focus2}.")


class Wedge(PhitsObject):
    name = "surface"
    syntax = surface_common | {"tip": ((None, None, None), (Real(), Real(), Real()), 0),
                       "sides": (None, OrthogonalMatrix(), 1)}


    shape = lambda self: ((f"*{self.index}" if self.reflective else
                           (f"+{self.index}" if self.white else f"{self.index}"),
                           "transform", "WED", "tip", " ".join(" ".join(str(i) for i in j) for j in self.sides)),)



class TetrahedronBox(PhitsObject):
    name = "surface"
    syntax = surface_common | {"xrange": ((None, None), Interval(), 0),
                               "yrange": ((None, None), Interval(), 1),
                               "zrange": ((None, None), Interval(), 2)}

    shape = lambda self: ((f"*{self.index}" if self.reflective else
                            (f"+{self.index}" if self.white else f"{self.index}"),
                            "transform", "RPP", "xrange", "yrange", "zrange"),)



_surface_spec = OneOf(IsA(Plane, index=True), IsA(PointPlane, index=True), IsA(ParallelPlane, index=True),
                     IsA(Sphere, index=True), IsA(Cylinder, index=True), IsA(Cone, index=True), IsA(SimpleConic, index=True),
                     IsA(GeneralConic, index=True), IsA(Box, index=True), # IsA(Torus, index=True), IsA(HexagonalPrism, index=True),
                     IsA(EllipticalCylinder, index=True), IsA(Spheroid, index=True), IsA(Wedge, index=True), IsA(TetrahedronBox, index=True))

## MISC

class MagneticField(PhitsObject):
    """A uniform magnetic field in a region, affecting charged particles."""
    name = "magnetic_field"
    syntax = {"typ": (None, FinBij({"dipole": 2, "quadrupole": 4}), 0),
              "strength": (None, Real(), 1),
              "calc_freq": (None, PosReal(), None, 0.0),
              "transform": (None, IsA(Transform, index=True), None, 0),
              "time": (None, PosReal(), None, "non"),
              }
    superobjects = ["cell"]
    prelude = (("reg", "'typ", "'gap", "mgf", "trcl", "'time"),)
    shape = (("cell", "typ", "calc_freq", "strength", "transform", "time"),)

    group_by = lambda self: type(self).__name__
    separator = lambda self: self.section_title()

class NeutronMagneticField(PhitsObject):
    """A uniform magnetic field in a region, affecting neutrons via spin."""
    name = "magnetic_field"
    syntax = {"typ": (None, FinBij({"identified": 60, "nograv": 61, "dipole": 102,
                                    "quadrupole": 104, "sextupole": 106}), 0),
              "strength": (None, Real(), 1),
              "calc_freq": (None, PosReal(), None, 0.0),
              "polarization": (None, Real(), None, "non"),
              "transform": (None, IsA(Transform, index=True), None, 0),
              "time": (None, PosReal(), None, "non"),
              }
    superobjects = ["cell"]
    prelude = (("reg", "'typ", "'gap", "mgf", "trcl", "polar", "'time"),)
    shape = (("cell", "typ", "calc_freq", "strength", "transform", "polarization", "time"),)

    group_by = lambda self: type(self).__name__
    separator = lambda self: self.section_title()


class MappedMagneticField(PhitsObject):
    """A non-uniform magnetic field in a region, affecting charged particles, given by a file."""
    name = "magnetic_field"
    syntax = {"typ": (None, FinBij({"xyz_list_charged": -1, "rz_list_charged": -2, "xyz_map_charged": -3, "rz_map_charged": -4,
                                    "xyz_list_neutron": -101, "rz_list_neutron": -102, "xyz_map_neutron": -103, "rz_map_neutron": -104}), 0),
              "normalization": (None, Real(), 1),
              "calc_freq": (None, PosReal(), 2),
              "m_file": (None, Path(), 3),
              "transform": (None, IsA(Transform, index=True), None, 0),
              }
    superobjects = ["cell"]
    prelude = (("reg", "'typ", "gap", "mgf", "trcl", "file"),)
    shape = (("cell", "typ", "calc_freq", "normalization", "transform", "m_file"),)

    group_by = lambda self: type(self).__name__
    separator = lambda self: self.section_title()



class ElectromagneticField(PhitsObject):
    """A uniform electromagnetic field within a region, affecting charged particles."""
    name = "electromagnetic_field"
    syntax = {"e_strength": (None, Real(), 0),
              "m_strength": (None, Real(), 1),
              "e_transform": (None, IsA(Transform, index=True), None, 0),
              "m_transform": (None, IsA(Transform, index=True), None, 0),
              }
    superobjects = ["cell"]
    prelude = (("reg", "elf", "mgf", "trcle", "trclm"),)
    shape = (("cell", "e_strength", "m_strength", "e_transform", "m_transform"),)

    group_by = lambda self: type(self).__name__
    separator = lambda self: self.section_title()



class MappedElectromagneticField(PhitsObject):
    """A non-unifrom electromagnetic field within a region, affecting charged particles, given by a file."""
    name = "electromagnetic_field"
    syntax = {"typ_e": (None, FinBij({"xyz_list_charged": -1, "rz_list_charged": -2, "xyz_map_charged": -3, "rz_map_charged": -4}), 0),
              "typ_m": (None, FinBij({"xyz_list_charged": -1, "rz_list_charged": -2, "xyz_map_charged": -3, "rz_map_charged": -4}), 1),
              "calc_freq": (None, PosReal(), 2),
              "e_normalization": (None, Real(), 3),
              "m_normalization": (None, Real(), 4),
              "e_file": (None, Path(), 5),
              "m_file": (None, Path(), 6),
              "e_transform": (None, IsA(Transform, index=True), None, 0),
              "m_transform": (None, IsA(Transform, index=True), None, 0)}
    superobjects = ["cell"]
    prelude = (("reg", "type", "typm", "gap", "elf", "mgf", "trcle", "trclm", "filee", "filem"),)
    shape = (("cell", "typ_e", "typ_m", "calc_freq", "e_normalization", "m_normalization", "e_transform", "m_transform", "e_file", "m_file"),)

    group_by = lambda self: type(self).__name__
    separator = lambda self: self.section_title()



class DeltaRay(PhitsObject):
    """A threshold energy fo each cell, above which delta rays are to be produced."""
    name = "delta_ray"
    syntax = {"threshold": (None, RealBetween(1, None), 0),
              }
    superobjects = ["cell"]
    prelude = (("reg", "del"),)
    shape = (("cell", "threshold"),)



class TrackStructure(PhitsObject):
    """Cell-by-cell setting of the track-structure model used for electrons/positrons."""
    name = "track_structure"
    syntax = {"model": (None, FinBij({"none": 0, "general": -1, "optimized": 1}), 0)}
    superobjects = ["cell"]
    prelude = (("reg", "mID"),)
    shape = lambda self: (("cell", "model"),)




class ElasticOption(PhitsObject):
    """Parameters for the user-specified elastic scattering law Fortran subroutines."""
    name = "elastic_option"
    syntax = {"c1": (None, PosReal(), 1),
              "c2": (None, PosReal(), 2),
              "c3": (None, PosReal(), 3),
              "c4": (None, PosReal(), 4)}

    prelude = (("reg", "'c1", "'c2", "'c3", "'c4"),)
    shape = (("cell", "c1", "c2", "c3", "c4"),)
    superobjects = ["cell"]



class FragData(PhitsObject):
    """Enables user-defined cross-sections for a particular interaction in a cell."""
    name = "frag_data"
    syntax = {"semantics": (None, FinBij({"histogram": 1, "extrapolated": 4, "interpolated": 5}), 0),
              "projectile": (None, OneOf(FinBij({"proton": "proton", "neutron": "neutron"}), Nuclide(fake=True)), 1),
              "target": (None, Nuclide(fake=True), 2),
              "file": (None, Path(), 3)}
    superobjects = ["cell"]
    prelude = (("opt", "proj", "targ", "'file"),)
    shape = (("semantics", "projectile", "target", "file"),)



class Importance(PhitsObject):
    """Change the relative tally weight of certain particles in a cell."""
    name = "importance"
    syntax = {"particles": ("part", List(Particle(), unique=True), 0),
              "importance": (None, PosReal(), 1),
              }
    superobjects = ["cell"]
    prelude = ("particles", ("reg", "imp"))
    shape = (("cell", "importance"),)
    group_by = lambda self: self.particles
    separator = lambda self: self.section_title()
    max_groups = 6

    @classmethod
    def global_restrictions(self, type_divided):
        all_particles = list(it.chain.from_iterable(map(lambda x: x.particles, type_divided["importance"])))
        if len(set(all_particles)) < len(all_particles):
            raise ValueError("Integration problem: all Importances must have mutually disjoint lists of particles.")


class WeightWindow(PhitsObject):
    """Makes the tally weight of some particle(s) in a region a function of time or energy."""
    name = "weight_window"
    syntax = {"particles": ("part", List(Particle(), unique=True), 0),
              "variable": (None, FinBij({"energy": "energy", "time": "time"}), 2),
              "windows": (None, List(Tuple(PosReal(), PosReal())), 1),
              }
    superobjects = ["cell"]
    prelude = lambda self: ("mesh = reg", "particles",
                            f"eng = {len(self.windows)}" if self.variable == "energy" else f"tim = {len(self.windows)}",
                            " ".join(map(lambda t: str(t[0]), self.windows)),
                            ("reg", " ".join(f"ww{i}" for i in range(1, len(self.windows) + 1))))
    shape = lambda self: (("cell", " ".join(map(lambda t: str(t[1]), self.windows))),)
    group_by = lambda self: (self.particles, self.variable)
    separator = lambda self: self.section_title()
    max_groups = 6

    @classmethod
    def global_restrictions(self, type_divided):
        all_particles = list(it.chain.from_iterable(map(lambda x: x.particles, type_divided["weight_window"])))
        if len(set(all_particles)) < len(all_particles):
            raise ValueError("Integration problem: all WeightWindows must have mutually disjoint lists of particles.")



class WWBias(PhitsObject):
    """Some magic with regards to the "variance reduction tecnique."""
    name = "ww_bias"
    syntax = {"particles": ("part", List(Particle(), unique=True), 0),
              "biases": (None, List(Tuple(PosReal(), PosReal())), 1),
              }
    superobjects = ["cell"]
    prelude = lambda self: ("particles", f"eng = {len(self.biases)}", " ".join(map(lambda t: str(t[0]), self.biases)),
                            ("reg", " ".join(f"wwb{i}" for i in range(1, len(self.biases) + 1))))
    shape = lambda self: (("cell", " ".join(map(lambda t: str(t[1]), self.biases))),)
    group_by = lambda self: self.particles
    separator = lambda self: self.section_title()
    max_groups = 6






class ForcedCollisions(PhitsObject):
    """Changes the way tallies are calculated in a region for better measurement of low-probability interactions, \
    such as thin targets, or for improving statistics."""
    name = "forced_collisions"
    syntax = {"particles": ("part", List(Particle(), unique=True), 0),
              "factor": (None, RealBetween(-1, 1), 1),
              "force_secondaries": (None, FinBij({True: 1, False: -1}), None),
              }

    superobjects = ["cell"]
    prelude = ("particles", ("reg", "fcl"))
    shape = lambda self: (("cell", f"{self.force_secondaries * self.factor}" if self.force_secondaries is not None \
                           else str(self.factor)),)

    group_by = lambda self: self.particles

    separator = lambda self: self.section_title()
    max_groups = 6

    def restrictions(self):
        if "electron" in self.particles or "positron" in self.particles:
            raise ValueError(f"ForcedCollision does not accept electrons or positrons; got {self.particles}.")

    def global_restrictions(self, type_divided):
        all_particles = list(it.chain.from_iterable(map(lambda x: x.particles, type_divided["forced_collisions"])))
        if len(set(all_particles)) < len(all_particles):
            raise ValueError("Integration problem: all ForcedCollisions must have mutually disjoint lists of particles.")



class RepeatedCollisions(PhitsObject):
    """Similar to `ForcedCollisions`, changes tally calculation in a region for low-probability interactions, \
    but with an eye towards rare, secondary-particle-producing reactions."""
    name = "repeated_collisions"
    syntax = {"particles": ("part", List(Particle(fake=True), unique=True), 0),
              "collision_reps": (None, PosInt(), 1),
              "evaporation_reps":  (None, PosInt(), 2),
              "mother": (None, List(Nuclide(fake=True)), 3),
              "ebounds": (("emin", "emax"), (PosReal(), PosReal()), None),

              }

    superobjects = ["cell"]
    prelude = lambda self: ("particles",
                            f"mother = {len(self.mother)}" if self.mother else "",
                            ("mother",),
                            "ebounds", ("reg", "n-coll", "n-evap"))
    shape = (("cell", "collision_reps", "evaporation_reps"),)

    group_by = lambda self: (self.particles, self.mother)
    separator = lambda self: self.section_title()
    max_groups = 6

    def restrictions(self):
        if self.collision_reps * self.evaporation_reps <= 1 or self.collision_reps * self.evaporation_reps >= 2_147_483_647:
            raise ValueError(f"RepeatedCollisions' product of repititions must be more than 1 as an int32;"
                             f" got collsion_reps={self.collision_reps} and evaporation_reps={self.evaporation_reps}.")
        if self.ebounds is not None and self.ebounds[0] >= self.ebounds[1]:

            raise ValueError(f"RepeatedCollisions' ebounds must be ordered; got {self.ebounds}.")

    @classmethod
    def global_restrictions(self, type_divided):
        for rc in type_divided["repeated_collisions"]:
            possible = set(map(lambda x: kf_encode(x[0]), rc.cell.material.composition))
            if any(kf_encode(x) not in possible for x in rc.mother):
                raise ValueError(f"Integration problem: RepeatedCollisions' mother nuclei must be among its cell's material's nuclei;"
                                 f" got {set(rc.mother) - possible} extra.")


# Relevant only for tallies we don't support.
# class Multiplier(PhitsObject):
#     name = "multiplier"
#     syntax = {"particles": ("part", List(Particle(), unique=True, max_len=6), 0),
#               "semantics": ("interpolation", FinBij({"linear": "lin", "log": "log", "left_histogram": "glow",
#                                                      "right_histogram": "ghigh"}), 1),
#               "bins": (None, List(Tuple(PosReal(), PosReal())), 2)}
#     shape = lambda self: (f"number = -{200 + self.index}", "semantics", "particles", f"ne = {len(self.bins)}",
#                           "\n".join(map(lambda t: f"{t[0]} {t[1]}", self.bins)))



class RegionName(PhitsObject):
    """Names a region in graphical output. Useful when `make_input`ing to visualise geometries."""
    name = "reg_name"
    syntax = {"reg_name": (None, Text(), 0),
              "size": (None, PosReal(), 1),
              }
    superobjects = ["cell"]
    shape = (("cell", "reg_name", "size"),)



# TODO: optional arguments?
class Counter(PhitsObject):
    """Configures one of three counters, which can track all manner of things. Results accessible through tallies."""
    name = "counter"
    syntax = {"particles": ("part", List(OneOf(Particle(), Nuclide()), max_len=20, unique=True), 0),
              "entry": (None, Between(-9999, 10000), 1),
              "exit": (None, Between(-9999, 10000), 2),
              "collision": (None, Between(-9999, 10000), 3),
              "reflection": (None, Between(-9999, 10000), 4),
              }

    superobjects = ["cell"]
    prelude = ("particles", ("reg", "in", "out", "coll", "ref"))
    shape = (("cell", "entry", "exit", "collision", "reflection"),)

    group_by = lambda self: self.particles
    separator = lambda self: self.section_title() + f"counter = {self.index}\n"
    max_groups = 3



class Timer(PhitsObject):
    """Configures the way time-of-flight is calculated within a region"""
    name = "timer"
    syntax = {"entry": (None, FinBij({"zero": -1, "nothing": 0, "stop": 1}), 1),
              "exit": (None, FinBij({"zero": -1, "nothing": 0, "stop": 1}), 2),
              "collision": (None, FinBij({"zero": -1, "nothing": 0, "stop": 1}), 3),
              "reflection": (None, FinBij({"zero": -1, "nothing": 0, "stop": 1}), 4),
              }
    superobjects = ["cell"]
    prelude = (("reg", "in", "out", "coll", "ref"),)
    shape = (("cell", "entry", "exit", "collision", "reflection"),)

## MATERIAL
_tester = Nuclide()
def _decomposition(composition):
    """Turns a composition list into PHITS input form."""
    r = ""
    for nuc, ratio in composition:
        conv = _tester.phits(nuc)
        if callable(conv):
            raise conv("composition")
        else:
            r += f"{conv} {ratio} "
    return r

class DataMax(PhitsObject):
    """Sets the maximum energy for an interaction between particles and nucleus in the material."""
    name = "data_max"
    syntax = {"particles": ("part", List(FinBij({"proton": "proton", "neutron": "neutron"}), unique=True), 0),
              "nucleus": (None, Nuclide(fake=True), 1),
              "max_energy": (None, PosReal(), 2)}
    superobjects = ["material"]
    prelude = ("particles", ("mat", "'nucleus", "dmax"))
    # Manual lies about accepting the usual syntax for nuclides
    shape = lambda self: (("material", "nucleus", "max_energy"),)
    group_by = lambda self: (self.particles,)
    separator = lambda self: self.section_title()
    max_groups = 6

class MatNameColor(PhitsObject):
    """Sets a name and color of a material. Useful when `make_input`ing to visualise geometries."""
    name = "mat_name_color"
    syntax = {"mat_name": (None, Text(), 0),
              "size": (None, PosReal(), 1),
              "color": (None, Color(), 2)} # TODO: color

    superobjects = ["material"]
    prelude = (("mat", "name", "'size", "'color"),)
    shape = (("material", "mat_name", "size", "color"),)


# TODO: temporarily setting the composition to be integer-only to cut down on length.
# TODO: also temporarily setting JENDL4Nuclide() to enable testing; it looks like Parameters() can change this error behavior,
# so should be set back to Nuclide() in production
# TODO: Nuclide-by-nuclide library setting
class Material(PhitsObject): # Composition is a list of pairs of (<element name string>, <ratio>) e.g. ("8Li", 0.5)
    """Sets up a material with which a `Cell` will be filled, including all relevant nuclear information."""
    name = "material"
    syntax = {"composition": (None, List(Tuple(JENDL4Nuclide(), PosInt()), unique_by=lambda x: kf_encode(x[0])), 0),
              "gas": ("GAS", Choice10(), None),
              "electron_step": ("ESTEP", PosInt(), None), # TODO: check integer right
              "neutron_lib": ("NLIB", LibraryID(), None),
              "photon_lib": ("PLIB", LibraryID(), None),
              "electron_lib": ("ELIB", LibraryID(), None),
              "proton_lib": ("HLIB", LibraryID(), None),
              "conductor": ("COND", FinBij({False: -1, True: 1}), None),
              "thermal_lib": (None, ThermalLib(), None),
              "chemical": ("chem", List(Tuple(Chemical(), PosInt())), None),
              "data_max": (None, IsA(DataMax, index=True), None),
              "mat_name_color": (None, IsA(MatNameColor, index=True), None)
              }
    subobjects = ["data_max", "mat_time_change", "mat_name_color"]
    shape = lambda self: (f"MAT[{self.index}]",
                          _decomposition(self.composition),
                          "gas", "electron_step", # "neutron_lib", "photon_lib", "electron_lib", "proton_lib",
                          "conductor",
                          "chem = " + " ".join(ch + " " + str(den) for ch, den in self.chemical) \
                          if self.chemical is not None else "",
                          f"MT{self.index} " + self.syntax["thermal_lib"][1].phits(self.thermal_lib) \
                          if self.thermal_lib is not None else ""
                          )

    def restrictions(self):
        if any(map(lambda x: int(str(_tester.phits(x[0]))[:-3]) > 92, self.composition)) and \
           (not hasattr(self, "parameters") or not hasattr(self.parameters, "stoping_model") or \
            "ATIMA" in self.parameters.stopping_model):
            raise ValueError("Material cannot have nuclei with Z > 92 and ATIMA set at the same time;"
                             " please pass stopping_model=SPAR+NTMC to the material.")


## CELL

_subobject_syntax = {"magnetic_field": (None, OneOf(IsA(MagneticField, index=True), IsA(NeutronMagneticField, index=True),
                                                   IsA(MappedMagneticField, index=True)), None),
                    "electromagnetic_field": (None, OneOf(IsA(ElectromagneticField, index=True),
                                                          IsA(MappedElectromagneticField, index=True)), None),
                    "delta_ray": (None, IsA(DeltaRay, index=True), None),
                    "track_structure": (None, IsA(TrackStructure, index=True), None),
                    "elastic_option": (None, IsA(ElasticOption, index=True), None),
                    "importance": (None, IsA(Importance, index=True), None),
                    "weight_window": (None, IsA(WeightWindow, index=True), None),
                    "ww_bias": (None, IsA(WWBias, index=True), None),
                    "forced_collisions": (None, IsA(ForcedCollisions, index=True), None),
                    "repeated_collisions": (None, IsA(RepeatedCollisions, index=True), None),
                    "reg_name": (None, IsA(RegionName, index=True), None),
                    "counter": (None, IsA(Counter, index=True), None),
                    "timer": (None, IsA(Timer, index=True), None)}

_cell_common_syntax = _subobject_syntax | {"volume": ("VOL", PosReal(), None),
                                           "temperature": ("TMP", PosReal(), None),
                                           "transform": ("TRCL", IsA(Transform, index=True), None)}

class Tetrahedral(PhitsObject):
    """A `Cell` that's a box filled with tetrahedrons from a file. Extremely computationally efficient."""
    name = "cell"
    syntax = _cell_common_syntax | {"regions": (None, IsA(TetrahedronBox, index=True), 0),
                              "material": (None, IsA(Material, index=True), 1),
                              "density": (None, PosReal(), 2),
                              "tet_format": (None, FinBij({"tetgen": "tetgen", "NASTRAN": "NASTRAN"}), 1),
                              "tet_file": (None, Path(), 2),
                              "scale_factor": ("TSFAC", PosReal(), None)}

    shape = lambda self: (("self", "material", "density", "regions", "\\"),
                          "volume\\", "temperature\\", "transform\\", "LAT=3\\",
                          f"tfile={self.tet_file}" if self.tet_format == "tetgen" else f"nfile={self.tet_file}", "scale_factor")

    subobjects = set(_subobject_syntax.keys())

    # def restrictions(self):
    #     if len(self.regions) != 1:
    #         raise ValueError(f"Tetrahedral cells may have only one TetrahedronBox region; got {self.regions}")
        # if self.forced_collisions is not None and self.repeated_collisions is not None:
        #     raise ValueError(f"Cannot set both forced_collisions and repeated_collisions on a Tetrahedral cell.")

    def __or__(self, other): # Union of cells; adopts leftmost's properties
        r = deepcopy(self)
        setattr(r, "regions", (self.regions,) + ("|",) + (other.regions,))
        return r

    def __invert__(self): # Set complement of cell; new cell has old properties
        r = deepcopy(self)
        r.regions = ("~", (self.regions,))
        return r

    def __and__(self, other): # Intersection of cells; drops properties
        r = deepcopy(self)
        r.regions = (self.regions,) + (other.regions,)
        return r
    def __rshift__(self, other): # returns other's regions with self's properties
        r = deepcopy(self)
        r.regions = other.regions
        return r

    def __lshift__(self, other): # returns self's region with other's properties
        r = deepcopy(other)
        r.regions = self.regions
        return r



class Void(PhitsObject):
    """A `Cell` with no material, just vacuum."""
    name = "cell"
    syntax = _cell_common_syntax | {"regions": (None, RegionTuple(_surface_spec), 0)}
    shape = lambda self: (("self", "0", "regions", "\\"), "volume\\", "temperature\\", "transform\\", "")
    subobjects = set(_subobject_syntax.keys())

    def __or__(self, other): # Union of cells; adopts leftmost's properties
        r = deepcopy(self)
        setattr(r, "regions", (self.regions,) + ("|",) + (other.regions,))
        return r

    def __invert__(self): # Set complement of cell; new cell has old properties
        r = deepcopy(self)
        r.regions = ("~", (self.regions,))
        return r

    def __and__(self, other): # Intersection of cells; drops properties
        r = deepcopy(self)
        r.regions = (self.regions,) + (other.regions,)
        return r

    def __rshift__(self, other): # returns other's regions with self's properties
        r = deepcopy(self)
        r.regions = other.regions
        return r

    def __lshift__(self, other): # returns self's region with other's properties
        r = deepcopy(other)
        r.regions = self.regions
        return r



class OuterVoid(PhitsObject):
    """Void, but different for some reason. Probably shouldn't be used directly;
    `run_phits` creates the required OuterVoid for you automatically."""
    name = "cell"
    syntax = _cell_common_syntax | {"regions": (None, RegionTuple(_surface_spec), 0)}
    shape = lambda self: (("self", "-1", "regions", "\\"), "volume\\", "temperature\\", "transform\\", "")
    subobjects = set(_subobject_syntax.keys())

    def __or__(self, other): # Union of cells; adopts leftmost's properties
        r = deepcopy(self)
        setattr(r, "regions", (self.regions,) + ("|",) + (other.regions,))
        return r

    def __invert__(self): # Set complement of cell; new cell has old properties
        r = deepcopy(self)
        r.regions = ("~", (self.regions,))
        return r

    def __and__(self, other): # Intersection of cells; drops properties
        r = deepcopy(self)
        r.regions = (self.regions,) + (other.regions,)
        return r
    def __rshift__(self, other): # returns other's regions with self's properties
        r = deepcopy(self)
        r.regions = other.regions
        return r

    def __lshift__(self, other): # returns self's region with other's properties
        r = deepcopy(other)
        r.regions = self.regions
        return r


class Cell(PhitsObject):
    """The prototypical `Cell`, consisting of the intersection of several regions defined by surfaces with a material of some density."""
    name = "cell"
    syntax = _cell_common_syntax | {"regions": (None, RegionTuple(_surface_spec), 0),
                              "material": (None, IsA(Material, index=True), 1),
                              "density": (None, PosReal(), 2)}
    shape = lambda self: (("self", "material", "density", "regions", "\\"),
                          "volume\\", "temperature\\", "transform\\", "")

    subobjects = set(_subobject_syntax.keys())

    def __or__(self, other): # Union of cells; adopts leftmost's properties
        r = deepcopy(self)
        setattr(r, "regions", (self.regions,) + ("|",) + (other.regions,))
        return r

    def __invert__(self): # Set complement of cell; new cell has old properties
        r = deepcopy(self)
        r.regions = ("~", (self.regions,))
        return r

    def __and__(self, other): # Intersection of cells; drops properties
        r = deepcopy(self)
        r.regions = (self.regions,) + (other.regions,)
        return r

    def __rshift__(self, other): # returns other's regions with self's properties
        r = deepcopy(self)
        r.regions = other.regions
        return r

    def __lshift__(self, other): # returns self's region with other's properties
        r = deepcopy(other)
        r.regions = self.regions
        return r

_cell_spec = OneOf(IsA(Cell, index=True), IsA(Tetrahedral, index=True), IsA(Void, index=True))
# idea: generate a UUID for the universe/fill, and then map UUIDs -> index at runtime
# other idea: make a Universe class, define an __init__, and make a call to super() for the normal __init__,
# but use the rest of __init__ to set the right attributes on the underlying cells, and marshall definitions
# def fill_universe(mask: Cell, contents: list[Cell]):
#     pass

## MISC 2

class SuperMirror(PhitsObject):
    """Enables calculation of low-energy neutron super-mirror reflections off the boundary between two `Cell`s via an empirical formula.\
    """
    name = "super_mirror"
    syntax = {"into": ("r-in", IsA(Cell, index=True), 0),
              "from": ("r-out", IsA(Cell, index=True), 1),
              "reflection_surface": ((None, None), (_surface_spec, _surface_spec), 2),
              "material_constant": (None, Real(), 3),
              "reflectivity": (None, Real(), 4),
              "critical_q": (None, Real(), 5),
              "falloff_rate": (None, Real(), 6),
              "cutoff_width": (None, PosReal(), 7)}

    prelude = (("r-in", "r-out", "mm", "r0", "qc", "am", "wm"),)
    shape = (("reflection_surface", "material_constant", "reflectivity", "critical_q", "falloff_rate", "cutoff_width"),)

class MatTimeChange(PhitsObject):
    """At a certain time, change the material old to material new."""
    name = "mat_time_change"
    syntax = {"time": (None, PosReal(), 0),
              "new": (None, IsA(Material, index=True), 1),
              "old": (None, IsA(Material, index=True), 2)}
    prelude = (("mat", "'time", "change"),)
    shape = (("old", "time", "new"),)



## SOURCE

# TODO: global scaling factor totfact, and correlation option iscorr. Something with group_by?

# removed from projectile spec: FinBij({"all": "all"})
_source_common = {"projectile": ("proj", List(OneOf(Particle(), Nuclide())), 0),
          "spin": (("sx", "sy", "sz"), (PosReal(), PosReal(), PosReal()), None),
          "mask": (("reg", "ntmax"), (IsA(Cell, index=True), PosInt()), None),
          "transform": ("trcl", IsA(Transform, index=True), None),
          "weight": ("wgt", PosReal(), None),
          "charge_override": ("izst", PosReal(), None),
          "counter_start": (("cnt(1)", "cnt(2)", "cnt(3)"), (PosInt(), PosInt(), PosInt()), None),
          "fissile": ("ispfs", FinBij({False: 0, "fissions": 1, "neutrons": 2}), None)
          # ibatch?
          }

_source_semi_common = {"elevation": ("dir", OneOf(RealBetween(0.0, 1.0), FinBij({"isotropic": "all"}), IsA(AngleDistribution)), None),
               "azimuth": ("phi", PosReal(), None),
               "dispersion": ("dom", OneOf(PosReal(), FinBij({"cos^2": -1})), None),
               # "energy": ("e0", PosReal(), 1), unsupported; just use a uniform energy distribution
               "spectrum": (None, IsA(EnergyDistribution), 1)}



class Cylindrical(PhitsObject):
    """A cylindrical solid source."""
    name = "source"
    syntax = _source_common | {"center": (("x0", "y0"), (Real(), Real()), None),
                        "zbounds": (("z0", "z1"), (Real(), Real()), None),
                        "radius": ("r0", PosReal(), None),
                        "cutout_radius": ("r1", PosReal(), None)} | _source_semi_common

    shape = lambda self: ("s-type = 1", "projectile", "spin", "mask", "transform", "weight", "charge_override", "counter_start",
                          "fissile", "center", "zbounds", "radius", "cutout_radius",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))

    def restrictions(self):
        if (self.radius is None or self.radius == 0) and self.cutout_radius is not None:
            raise ValueError("Cylindrical sources that specify a cutout radius must also specify a nonzero radius;"
                             "got cutout_radius={self.cutout_radius}.")
        if self.radius is not None and self.cutout_radius is not None and self.radius < self.cutout_radius:
            raise ValueError("Cylindrical sources cannot have cutouts larger than their radius;"
                             f"got radius={self.radius} and cutout_radius={self.cutout_radius}.")

class Rectangular(PhitsObject):
    """A rectangular solid source."""
    name = "source"
    syntax = _source_common | {"xbounds": (("x0", "x1"), (Real(), Real()), None),
                       "ybounds": (("x0", "x1"), (Real(), Real()), None),
                       "zbounds": (("x0", "x1"), (Real(), Real()), None)} | _source_semi_common

    shape = lambda self: ("s-type = 2", "projectile", "spin", "mask", "transform", "weight", "charge_override", "counter_start",
                          "fissile", "xbounds", "ybounds", "zbounds",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))




class Gaussian(PhitsObject):
    """A Gaussian source from every direction."""
    name = "source"
    syntax = _source_common | {"center": (("x0", "y0", "z0"), (Real(), Real(), Real()), None),
                       "fwhms": (("x1", "y1", "z1"), (PosReal(), PosReal(), PosReal()), None)} | _source_semi_common

    shape = lambda self: ("s-type = 3", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "center", "fwhms",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                          else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))

class GaussianSlices(PhitsObject):
    """A 2D Gaussian source, uniform in the \\(z\\)-axis."""
    name = "source"
    syntax = _source_common | {"center": (("x0", "y0"), (Real(), Real()), None),
                       "fwhm": ("r1", PosReal(), None),
                       "zbounds": (("z0", "z1"), (Real(), Real()), None)} | _source_semi_common

    shape = lambda self: ("s-type = 13", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "center", "fwhm", "zbounds",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                          else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))




class Parabolic(PhitsObject):
    """A parabolic source from every direction."""
    name = "source"

    syntax = _source_common | {"center": (("x0", "y0"), (Real(), Real()), None),
                       "width": (("x1", "y1"), (PosReal(), PosReal()), None),
                       "zbounds": (("z0", "z1"), (Real(), Real()), None),
                       "order": ("rn", Between(2, 2147483647), None) # PHITS's default INTEGER is 32-bit; if something's bigger,
                                                                     # their ≡ 0 (mod 2) check of multiplying and dividing by 2 fails.
                       } | _source_semi_common
    shape = lambda self: ("s-type = 7", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "center", "width", "zbounds", "order",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))

    def restrictions(self):
        if self.order is not None and self.order % 2 != 0:
            raise ValueError(f"The order of a Parabolic source must be even; got order={self.order}.") # TODO: needed?
        if self.zbounds is not None and self.zbounds[0] > self.zbounds[1]:
            raise ValueError(f"The the zbounds of a Parabolic source must be a well-formed interval; got zbounds={self.zbounds}.")


class ParabolicSlices(PhitsObject):
    """A parabolic source, uniform along the \\(z\\)-axis."""
    name = "source"
    syntax = _source_common | {"center": (("x0", "y0"), (Real(), Real()), None),
                       "width": ("r1", Real(), None),
                       "zbounds": (("z0", "z1"), (Real(), Real()), None),
                       "order": ("rn", Between(2, 2147483647), None)
                       } | _source_semi_common
    shape = lambda self: ("s-type = 15", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "center", "width", "zbounds", "order",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))

    def restrictions(self):
        if self.order is not None and self.order % 2 != 0:
            raise ValueError(f"The order of a ParabolicPrism source must be even; got order={self.order}.") # TODO: needed?
        if self.zbounds is not None and self.zbounds[0] > self.zbounds[1]:
            raise ValueError(f"The the zbounds of a ParabolicPrism source must be a well-formed interval; got zbounds={self.zbounds}.")


# dir = iso not supported
class Spherical(PhitsObject):
    """A spherical or spherical-shell solid source."""
    name = "source"
    syntax = _source_common | {"center": (("x0", "y0", "z0"), (Real(), Real(), Real()), None),
                       "r_in": ("r1", PosReal(), None),
                       "r_out": ("r2", PosReal(), None),
                       # "elevation_bounds": (("ag1", "ag2"), (Real(), Real()), None),
                       # "azimuth_bounds": (("pg1", "pg2"), (Real(), Real()), None),
                       "elevation": ("dir", OneOf(RealBetween(0.0, 1.0), FinBij({"all": "all"}), IsA(AngleDistribution)), None),
                       # TODO: this elevation and this elevation only doesn't work if I set FinBij({"isotropic": "all"}).
                       "resample_cutoff": ("isbias", Choice10(), None),
                       "spectrum": (None, IsA(EnergyDistribution), 1)}
    shape = lambda self: ("s-type = 9", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "center", "r_in", "r_out",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "resample_cutoff", ("spectrum",))

    def restrictions(self):
        if (self.elevation == "isotropic" and self.r_in is not None and (self.r_out is None or self.r_out == 0)) \
           or (self.elevation == "isotropic" and self.r_in is not None and self.r_out is not None and self.r_in <= self.r_out):
            raise ValueError("Spherical sources with isotropic elevation must have greater inner radius than outer radius;"
                             f"got r_in={self.r_in} and r_out={self.r_out}.")


        if (self.elevation != "isotropic" and self.r_in is not None and (self.r_out is None or self.r_out == 0)) \
           or (self.elevation != "isotropic" and self.r_in is not None and self.r_out is not None and self.r_in > self.r_out):
            raise ValueError("Spherical sources that specify an inner radius must also specify a greater outer radius;"
                             f"got r_in={self.r_in} and r_out={self.r_out}.")


class Beam(PhitsObject): # I don't understand what this is trying to do
    """A beam-like source."""
    name = "source"
    syntax = _source_common | {"center": (("x0", "y0"), (Real(), Real()), None),
                       "eccentricity": (("x1", "y1"), (Real(), Real()), None),
                       "zbounds": (("z0", "z1"), (Real(), Real()), None),
                       "phase_gradients": (("rx", "ry"), (Real(), Real()), None),
                       "sampling": ("wem", OneOf(FinBij({"gaussian": 0}), PosReal()), None),
                       "dispersion": (("x1", "y1"), (Real(), Real()), None),
                       "angle_dispersion": (("xmrad1", "ymrad1"), (PosReal(), PosReal()), None),
                       "phase_center": (("x2", "y2"), (Real(), Real()), None),
                       "phase_angle_center": (("xmrad2", "ymrad2"), (Real(), Real()), None),
                       "positive": ("dir", FinBij({True: 1, False: -1}), None),
                       "spectrum": (None, IsA(EnergyDistribution), 1)}

    shape = ("s-type = 11", "projectile", "spin", "mask", "transform", "weight", "counter_start",
             "charge_override", "fissile", "center", "eccentricity", "zbounds", "phase_gradients", "sampling", "dispersion",
             "angle_dispersion", "phase_center", "phase_angle_center", "positive", ("spectrum",))


# decay-turtle??????


class Conical(PhitsObject):
    """A conical solid source."""
    name = "source"
    syntax = _source_common | {"top": (("x0", "y0", "z0"), (Real(), Real(), Real()), None),
                       "altitude": (("x1", "y1", "z1"), (Real(), Real(), Real()), None),
                       "trim": (("r0", "r1"), (Real(), Real()), None),
                       "angle": ("r2", PosReal(), None)} | _source_semi_common
    shape = lambda self: ("s-type = 18", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "top", "altitude", "trim", "angle",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth",
                          "dispersion", ("spectrum",))



class TriangularPrism(PhitsObject):
    """A triangular-prism solid source."""
    name = "source"
    syntax = _source_common | {"origin": (("x0", "y0", "z0"), (Real(), Real(), Real()), None),
                       "side1": (("x1", "y1", "z1"), (Real(), Real(), Real()), None),
                       "side2": (("x2", "y2", "z2"), (Real(), Real(), Real()), None),
                       "extrusion": (("x3", "y3", "z3"), (Real(), Real(), Real()), None),
                       "attenuation": ("exa", PosReal(), None)} | _source_semi_common
    shape = lambda self: ("s-type = 20", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "origin", "side1", "side2", "extrusion", "attenuation",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth",
                          "dispersion", ("spectrum",))


class TetrahedralSource(PhitsObject): # TODO: subobjects
    """A `Tetrahedral`ly-defined solid source."""
    name = "source"
    syntax = _source_common | {"cell": ("tetreg", IsA(Tetrahedral, index=True), 2)} | _source_semi_common
    shape = lambda self: ("s-type = 24", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "cell",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))


class SurfaceSource(PhitsObject):
    """A solid source defined by some part of a surface."""
    name = "source"
    syntax = _source_common | {"surface": ("suf", _surface_spec, 2),
                               "cut": ("cut", List(_surface_spec, max_len=8), 3)} | _source_semi_common
    shape = lambda self: ("s-type = 26", "projectile", "spin", "mask", "transform", "weight", "counter_start",
                          "charge_override", "fissile", "surface", "cut",
                          (f"dir = data\n{self.elevation.definition()}" if isinstance(self.elevation, AngleDistribution) \
                           else f"dir = {self.elevation}") if self.elevation is not None else "", "azimuth", "dispersion", ("spectrum",))

_source_spec = OneOf(IsA(Cylindrical, index=True), IsA(Rectangular, index=True), IsA(Gaussian, index=True),
                    IsA(GaussianSlices, index=True), IsA(Parabolic, index=True), IsA(ParabolicSlices, index=True),
                    IsA(Spherical, index=True), IsA(Beam, index=True), IsA(Conical, index=True), IsA(TriangularPrism, index=True),
                    IsA(TetrahedralSource, index=True), IsA(SurfaceSource, index=True))

# No dump file or user source (for the latter, you can write your own PhitsObject)
# TODO: dom = -10 for Cylindrical and Rectangular

## TALLY

# counters? multiplier?
class DumpFluence(PhitsObject):
    """Tally that counts some quantity whenever a particle crosses from one region to another."""
    name = "t-cross"
    syntax = {"out": (None, IsA(Cell, index=True), 0),
              "into": (None, IsA(Cell, index=True), 1),
              "area": (None, PosReal(), 2),
              "data": ("dump", List(FinBij({"particle": 1, "x": 2, "y": 3, "z": 4, "u": 5, "v": 6, "w": 7, "energy": 8, "weight": 9,
                                            "time": 10, "counter1": 11, "counter2": 12, "counter3": 13, "spinx": 14, "spiny": 15,
                                            "spinz": 16, "collision_number": 17, "history_number": 18, "batch_number": 19,
                                            "cascade_id": 20}), unique=True), 3),
              "output": ("output", FinBij({"current": "current", "a_current": "a-curr", "oa_current": "oa-curr"}), 4),
              "particles": ("part", List(Particle(), max_len=6, unique=True), None),
              "factor": ("factor", PosReal(), None),
              # "counter": ()
              "maximum_error": ("stdcut", PosReal(), None),
              # "multiplier": ()
              # TODO: set Between() back to PosInt; done because otherwise memory overflows on my machine during tests
              "energy_mesh": (("emin", "emax", "ne"), (Interval(0), Between(1, 50)), 5),
              "angle_mesh": (("amin", "amax", "na"), (Interval(-1, 1), Between(1, 50)), 6),
              "angle_semantics": ("iangform", FinBij({"to_normal": 0, "to_x": 1, "to_y": 2, "to_z": 3}), None),
              "time_mesh": (("tmin", "tmax", "nt"), (Interval(0), Between(1, 50)), None)}

    prelude = lambda self: ("particles", "unit = 1", "axis = reg", f"file = cross{self.index}", "factor", "output", "maximum_error",
                            f"dump = -{len(self.data)}", ("data",),
                            "e-type = 2", f"emin = {self.energy_mesh[0][0]}", f"emax = {self.energy_mesh[0][1]}",
                            f"ne = {self.energy_mesh[1]}",
                            f"a-type = 2\namin = {self.angle_mesh[0][0]}\namax = {self.angle_mesh[0][1]}\nna = {self.angle_mesh[1]}" \
                            if self.angle_mesh is not None else "", "angle_semantics",
                            f"t-type = 2\ntmin = {self.angle_mesh[0][0]}\ntmax = {self.time_mesh[0][1]}\nnt = {self.time_mesh[1]}" \
                            if self.time_mesh is not None else "",
                            "mesh = reg", f"reg = {self.group_size}",
                            ("r-from", "r-to", "'area"))

    shape = lambda self: ((f"{self.out.index}", f"{self.into.index}", "area"),)

    group_by = lambda self: (self.particles, self.data, self.output, self.factor, self.energy_mesh, self.angle_mesh,
                             self.time_mesh)
    separator = lambda self: self.section_title()


class DumpProduction(PhitsObject):
    """Tally that counts how many particles of some type are created within a `Cell`."""
    name = "t-product"
    syntax = {"cell": ("reg", IsA(Cell, index=True), 0),
              "data": ("dump", List(FinBij({"particle": 1, "x": 2, "y": 3, "z": 4, "u": 5, "v": 6, "w": 7, "energy": 8, "weight": 9,
                                            "time": 10, "counter1": 11, "counter2": 12, "counter3": 13, "spinx": 14, "spiny": 15,
                                            "spinz": 16, "collision_number": 17, "history_number": 18, "batch_number": 19,
                                            "cascade_id": 20}), unique=True), 1),
              "output": ("output", FinBij({"source": "source", "nuclear": "nuclear", "nonela": "nonela", "elastic": "elastic",
                                           "decay": "decay", "fission": "fission", "atomic": "atomic"}), 2),
              "particles": ("part", List(Particle(), max_len=6, unique=True), None),
              "materials": ("material", List(IsA(Material, index=True)), None),
              "mother": ("material", List(Nuclide(fake=True)), None),
              "factor": ("factor", PosReal(), None),
              # counter
              "maximum_error": ("stdcut", PosReal(), None),
              # "multiplier": ()
              # TODO: set Between() back to PosInt; done because otherwise the memory overflows
              "energy_mesh": (("emin", "emax", "ne"), (Interval(0), Between(1, 50)), 3),
              "time_mesh": (("tmin", "tmax", "nt"), (Interval(0), Between(1, 50)), None)}


    prelude = lambda self: ("particles", "unit = 1", "axis = reg", f"file = product{self.index}", "factor", "output", "maximum_error",
                            f"material = {len(self.materials)}" if self.materials is not None else "", ("materials",),
                            f"mother = {len(self.mother)}" if self.mother is not None else "", ("mother",),
                            f"dump = -{len(self.data)}", ("data",),
                            "e-type = 2", f"emin = {self.energy_mesh[0][0]}", f"emax = {self.energy_mesh[0][1]}",
                            f"ne = {self.energy_mesh[1]}",
                            f"t-type = 2\ntmin = {self.time_mesh[0][0]}\ntmax = {self.time_mesh[0][1]}\nnt = {self.time_mesh[1]}" \
                            if self.time_mesh is not None else "",
                            "mesh = reg")

    shape = ("cell",)

    group_by = lambda self: (self.particles, self.data, self.output, self.factor, self.energy_mesh, self.time_mesh)
    separator = lambda self: self.section_title()


class DumpTime(PhitsObject):
    """Tally that records how particles' disappearance (via energy cutoff, escape, decay) changes over time in a `Cell`."""
    name = "t-time"
    syntax = {"cell": ("reg", IsA(Cell, index=True), 0),
              "data": ("dump", List(FinBij({"particle": 1, "x": 2, "y": 3, "z": 4, "u": 5, "v": 6, "w": 7, "energy": 8, "weight": 9,
                                            "time": 10, "counter1": 11, "counter2": 12, "counter3": 13, "spinx": 14, "spiny": 15,
                                            "spinz": 16, "collision_number": 17, "history_number": 18, "batch_number": 19,
                                            "cascade_id": 20}), unique=True), 1),
              "output": ("output", FinBij({"all": "all", "cutoff": "cutoff", "escape": "escape", "decay": "decay"}), 2),
              "particles": ("part", List(Particle(), max_len=6, unique=True), None),
              "materials": ("material", List(IsA(Material, index=True)), None),
              "factor": ("factor", PosReal(), None),
              # counter
              "maximum_error": ("stdcut", PosReal(), None),
              # "multiplier": ()
              "energy_mesh": (("emin", "emax", "ne"), (Interval(0), Between(1, 50)), 3),
              "time_mesh": (("tmin", "tmax", "nt"), (Interval(0), Between(1, 50)), 4)}


    prelude = lambda self: ("particles", "unit = 1", "axis = reg", f"file = time{self.index}", "factor", "output", "maximum_error",
                            f"material = {len(self.materials)}" if self.materials is not None else "", ("materials",),
                            f"dump = -{len(self.data)}", ("data",),
                            "e-type = 2", f"emin = {self.energy_mesh[0][0]}", f"emax = {self.energy_mesh[0][1]}",
                            f"ne = {self.energy_mesh[1]}",
                            "t-type = 2", f"tmin = {self.time_mesh[0][0]}", f"tmax = {self.time_mesh[0][1]}",
                            f"nt = {self.time_mesh[1]}", "mesh = reg")

    shape = ("cell",)

    group_by = lambda self: (self.particles, self.data, self.output, self.factor, self.energy_mesh, self.time_mesh)
    separator = lambda self: self.section_title()

_tally_spec = OneOf(IsA(DumpFluence, index=True), IsA(DumpProduction, index=True), IsA(DumpTime, index=True))



## DMP-READER

def read_dump(name: str, columns: list[str], return_type: str) -> dict:
    """Given a path to a PHITS dump file and names of the record entries in order,
    produce a semantically equivalent dictionary of lists/numpy array/Pandas dataframe of the contents."""
    rd = FortranRecordReader('(30(1p1d24.15))') # PHITS documentation says e, but code says d---latter is consistent with behavior
    acc = dict.fromkeys(columns)
    with open(name, 'r') as dmp:
        for line in dmp:
            for i, val in enumerate(rd.read(line)):
                if val is not None:
                    col = columns[i]
                    if acc[col] is None:
                        acc[col] = [val if col != "particle" else kf_decode(val)]
                    else:
                        acc[col] += val if col != "particle" else kf_decode(val)

    if return_type == "dict":
        return acc
    elif return_type == "numpy":
        return np.fromiter(acc.items())
    elif return_type == "pandas":
        return DataFrame.from_dict(acc)


    return acc



## RUN-PHITS


def make_input(cells, sources, tallies, title: str = str(datetime.now()), cross_sections=[], multipliers=[], super_mirrors=[],
               mat_time_changes=[], raw="", outer_void_properties: dict = dict(), **kwargs) -> str:
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
    assert (isinstance(cells, PhitsObject) and cells.name == "cell") \
        or (len(cells) >= 1 and all(map(lambda x: x.name == "cell", cells))), \
    f"`cells` must be either a cell object or a list of cell objects; got {cells}."

    assert (isinstance(sources, PhitsObject) and sources.name == "source") \
        or (len(sources) >= 1 and all(map(lambda x: x.name == "source", sources))), \
    f"`sources` must be either a sources object or a list of source objects; got {sources}."

    assert (isinstance(tallies, PhitsObject) and "t-" in tallies.name) \
        or (all(map(lambda x: "t-" in x.name, tallies))), \
    f"`tallies` must be either a tally object or a list of tally objects; got {tallies}."

    assert (isinstance(cross_sections, PhitsObject) and cross_sections.name == "frag_data") \
        or (all(map(lambda x: x.name == "frag_data", cross_sections))), \
    f"`cross_sections` must be either a FragData object or a list of FragData objects; got {cross_sections}."

    assert (isinstance(super_mirrors, PhitsObject) and super_mirrors.name == "super_mirror") \
        or (all(map(lambda x: x.name == "super_mirror", super_mirrors))), \
    f"`super_mirror` must be either a SuperMirror object or a list of SuperMirror objects; got {super_mirrors}."

    assert (isinstance(mat_time_changes, PhitsObject) and mat_time_changes.name == "mat_time_change") \
        or (all(map(lambda x: x.name == "mat_time_change", mat_time_changes))), \
    "`mat_time_change` must be either a MatTimeChange object or a list of MatTimeChange objects; "
    f"got {mat_time_changes}."

    # Problem: you can have different objects that are "essentially the same" appearing in the object tree.
    # Solution: exploit the __eq__ and __hash__ defined on PhitsObject
    unique = set()

    def add_to_set(an_obj, the_set, prev=None):  # Recursively add subtypes to set if they represent an "entry" in one of the sections
        if isinstance(an_obj, list) or isinstance(an_obj, tuple):
            for ob in an_obj:
                if ob is not prev:
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
    add_to_set(super_mirrors, unique)
    add_to_set(mat_time_changes, unique)


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
                    "super_mirror": [],
                    "elastic_option": [],
                    "data_max": [],
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

    toset = OuterVoid((~reduce(lambda c1, c2: c1 | c2, type_divided["cell"])).regions, **outer_void_properties)

    type_divided["cell"].append(toset)


    for section, entries in type_divided.items():
        for idx, value in enumerate(entries):
            value.index = idx+1

    # Problem: while we've chosen a set of representatives for equivalence classes under PhitsObject.__eq__, the objects themselves
    # don't have subobjects with index attributes pointing to the representative---there will be None showing up all over the output.
    # Solution: replace all members of an equivalence class in the object tree with their representative (whose index is defined above).
    representatives = {n: n for n in it.chain.from_iterable(type_divided.values())} # necessary because `unique` doesn't have idx

    def replace(this, that, inside):
        if inside is this:
            return that
        elif isinstance(inside, list) or isinstance(inside, tuple):
            replaced = tuple()
            for el in inside:
                if el is this:
                    replaced += (that,)
                else:
                    replaced += (replace(this, that, el),)
            return replaced
        else:
            return this

    def adjust_subobjects(an_obj, ason=(None, None)): # Recursively replace redundant subtypes with the representative in the dict
        if isinstance(an_obj, tuple) or isinstance(an_obj, list):
            for ob in an_obj:
                if ob is not ason[1]:
                    adjust_subobjects(ob, ason=ason)

        elif isinstance(an_obj, PhitsObject):
            if ason != (None, None):
                representative = representatives[an_obj]
                to_check = getattr(ason[1], ason[0])
                setattr(ason[1], ason[0], replace(an_obj, representatives[an_obj], to_check))
                
            for name, child in an_obj.__dict__.items():
                if child is not ason[1]:
                    adjust_subobjects(child, ason=(name, an_obj))

    adjust_subobjects(cells)
    adjust_subobjects(sources)
    adjust_subobjects(tallies)
    adjust_subobjects(cross_sections)
    adjust_subobjects(multipliers)
    adjust_subobjects(super_mirrors)
    adjust_subobjects(mat_time_changes)


    # Check that the whole shebang is valid together
    for _, v in type_divided.items():
        if len(v) > 0 and hasattr(type(v[0]), "global_restrictions"):
            v[0].global_restrictions(type_divided)

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
    add_defs("super_mirror")
    add_defs("elastic_option")
    add_defs("data_max")
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

    if max(map(len, inp.split("\n"))) > 200:
        print(inp)
        raise RuntimeError("PHITS line limit reached.")
    else:
        return inp





def run_phits(cells, sources, tallies, command: str = "phits", hard_error: bool = True, filename: str = "phits.inp",
              return_type: str = "dict", injected_files=[], yanked_files=[], yank_to=None, **make_input_kwargs):
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
    | `injected_files` | A list of path- or file-like objects to `shutil.copy`/`shutil.copyfileobj` to the temporary directory, or a dict of filenames and contents to be inserted there. |
    | `yanked_files` | A list of paths relative to the temporary directory which will be extracted before cleanup. |
    | `yank_to` | A path to which `yanked_files` will be saved---if left `None`, the file contents will be returned as strings. |
    """
    with tf.TemporaryDirectory() as newdir:
        inp = make_input(cells, sources, tallies, **make_input_kwargs)
        name = os.path.join(newdir, filename)
        with open(os.path.join(newdir, filename), "w") as inp_file:
            inp_file.write(inp)

        if isinstance(injected_files, list): # names of files to copy
            for path in injected_files:
                copy(path, os.path.join(newdir, os.path.basename(path)))
        else: # dict of file names: contents to write
            assert isinstance(injected_files, dict), "injected_files must be either a list of filenames or a dict of filename: contents."
            for k, v in injected_files.items():
                with open(os.path.join(newdir, k), "w") as inj:
                    inj.write(v)

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
            result[t] = read_dump(dfile, t.data, return_type)


        if yanked_files != []:
            for name in yanked_files:
                if yank_to is not None:
                    copy(os.path.join(newdir, name), yank_to)
                else:
                    with open(os.path.join(newdir, name), "w") as contents:
                        result[name] = contents.read()

        return result

__pdoc__ = dict()
__pdoc__["builds"] = False
__pdoc__["slices"] = False
__pdoc__["builds_right"] = False
__pdoc__["interval"] = False
for name, cl in list(sys.modules[__name__].__dict__.items()):
    if type(cl) == type and issubclass(cl, PhitsObject) and cl != PhitsObject:
        __pdoc__[cl.__name__] = cl.__doc__ + cl.syntax_desc() if cl.__doc__ else cl.syntax_desc()

__pdoc__["Parameters"] = Parameters.__doc__ + Parameters.syntax_desc()
