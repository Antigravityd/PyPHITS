
from base import *


class MatTimeChange(PhitsObject):
    name = "mat_time_change"
    required = ["time", "new"]
    positional = ["time", "new"]
    optional = ["old"]
    shape = (("old", "time", "new"))
    prelude = (("mat", "\\time", "change"))
    nones = {"time": "non", "change": "non"}



class DataMax(PhitsObject): # requires special handling in make_input
    name = "data_max"
    required = ["particles", "nucleus", "threshold"]
    positional = ["particles", "nucleus", "threshold"]
    optional = ["material"]
    shape = (("material", "nucleus", "threshold"))
    prelude = ("particles", (("mat", "\\nucleus", "dmax")))
    group_by = lambda self: self.particles
    max_groups = 6
    separator = lambda self: self.section_title()
    ident_map = {"particles": "part"}


class MatNameColor(PhitsObject):
    name = "mat_name_color"
    required = ["name", "size", "color"]
    positional =  ["name", "size", "color"]
    optional = ["material"]
    shape = (("material", "name", "size", "color"))
    prelude = (("mat", "\\name", "\\size", "\\color"))



# TODO: how the libraries work isn't well-documented. Is there a single library set for the whole material, or
# does one set a library after each element of the compositon? Can the thermal neutron library be set anywhere?
# dedxfile?
# den nuc?
class Material(PhitsObject): # Composition is a list of pairs of (<element name string>, <ratio>) e.g. ("8Li", 0.5)
    name = "material"
    required = ["composition"]
    mapping = {"chemical": ("chem", None, ),
               "condensed": ("GAS", True, choice10(c_style=True)),
               "electron_substep": ("ESTEP", ),
               "neutron_lib": (),
               "proton_lib": (),
               "electron_lib": (),
               "photon_lib": (),
               "thermal_libs": ()}
    prefix_parser = r"""
    materials: "[" "m" "a" "t" "e" "r" "i" "a" "l" "]" "\n" columnorder?
    columnorder: "nuc" "den" | "den" "nuc"
    """
    parser = r"""
    %ignore WS

    start: ("MAT[" POSINT "]" | "M" POSINT) (comp_entry | assignment)+ [thermlib | chemical | thermlib chemical | chemical thermlib]


    comp_entry: element computation
    thermlib: "MT" POSINT (/[\.0-9a-z]+/)+
    chemical: "chem" "="  (CHEMSYM computation)+

    CHEMSYM: "H2O" | "CO2" | "NH2" | "NH3" | "SF6" | "TeF6" | "CH4" | "CH3" | "C2H2" | "C2H4" | "C2H6" | "C6H6" | "(CH3)2NH"
    element: (PREFIX_ELT | POSTFIX_ELT | KF_ELT) POSTFIX_LIB?
    assignment: IDENTIFIER "=" computation
    PREFIX_ELT: /[a-z]+-[0-9]+/
    POSTFIX_ELT: /[0-9]*[a-z]+/
    KF_ELT: /[0-9]+/
    POSTFIX_LIB: /\.[0-9][0-9][a-z]/"""

    POSTFIX_LIB: str
    KF_ELT: str
    POSTFIX_ELT: str
    PREFIX_ELT: str
    def assignment(self, ass):
        return Token("assignment", (ass[0], ass[1]))

    def element(self, elt):
        # TODO: we distinguish the postfix-element syntax e.g 208Pb. Consider.
        if elt[0].type == "PREFIX_ELT":
            split = str(elt[0]).split("-")
            return f"{split[1]}{split[0]}" if len(elt) == 1 else (f"{split[1]}{split[0]}", elt[1])
        elif elt[0].type = "POSTFIX_ELT":
            return str(elt[0]) if len(elt) == 1 else (str(elt[0]), elt[1])
        elif elt[0].type == "KF_ELT":
            # TODO: sym_lookup
            nice = f"{str(elt[0])[-3:]}{sym_lookup(str(elt[0])[-3:])}"
            return nice if len(elt) == 1 else (nice, elt[1])

    def chemical(self, chem):
        return Token("chemical", zip([str(i) for i in chem[::2]], chem[1::2]))

    def thermlib(self, lib):
        return Token("thermlib", [str(i) for i in lib[:1]])


    def start(self, mat):
        composition = [(str(el), float(wgt)) for el, wgt in filter(lambda tok: isinstance(tok, Token) and tok.type == "comp_entry", mat)]
        assignments = {self.ident_to_py(str(k)): float(v) for k, v in
                       filter(lambda tok: isinstance(tok, Token) and tok.type == "assignment", mat)}
        chemical = list(filter(lambda tok: isinstance(tok, Token) and tok.type == "chemical") for i in sublist)
        thermlibs = list(filter(lambda tok: isinstance(tok, Token) and tok.type == "thermlib"))

        return Material(composition, chemical=chemical, thermal_libs=thermlibs, **assignments)


    def tree(self):
        start = [Token("POSINT", self.index)]
        for elt, dens in self.composition:
            if len(elt) == 1:
                start.append(Tree("comp_entry", [Tree("element", [Token("POSTFIX_ELT", elt[0])]), Token("computation", dens)]))
            else:
                start.append(Tree("comp_entry", [Tree("element", [Token("POSTFIX_ELT", elt[0]), Token("POSTFIX_LIB", elt[1])]),
                                                 Token("computation", dens)]))
        for k, v in {k, self.__dict__[k] for k in ["condensed", "electron_substep", "neutron_lib", "electron_lib", "photon_lib"]}.items():
            start.append(Tree("assignment", [Token("IDENTIFIER", ident_to_phits(k)), Tree("computation", value_to_phits(k, v))]))



        start.append(Tree("thermlib", [Token("POSINT", self.index), self.thermal_libs]))
        start.append(Tree("chemical", [(Token("CHEMSYM", i[0]), Tree("computation", i[1])) for i in self.chemical]))

        return Tree("start", start)


    subobjects = ["time_change", "data_max", "mat_name_color"]


def chemsymbol()
