from base import *

class Rectangle(PhitsObject):
    name = "rectangle"
    syntax = {"width": ("tw", Real(), 0),
              "number": ("tn", PosInt(), 1),
              "delta": ("td", Real(), 2),
              "center": ("t0", Real(), None)}

    _grammar = r'''start:  "t-type" "=" "1" "\n" normal ~ 3..4
    normal: @assign_among|self.inv_syntax().keys()|
    %ignore SPACE
    '''

    def tree(self):
        tree = [Tree(Token("RULE", "normal"), [Tree(Token("RULE", "tw"), [Token("computation", self.width)])]),
                Tree(Token("RULE", "normal"), [Tree(Token("RULE", "tn"), [Token("POSINT", self.number)])]),
                Tree(Token("RULE", "normal"), [Tree(Token("RULE", "td"), [Token("computation", self.delta)])])]
        if self.center is not None:
            tree.append(Tree(Token("RULE", "normal"), [Tree(Token("RULE", "t0")), [Token("computation", self.center)]]))

        return Tree(Token("RULE", "start"), tree)



class Gaussian(PhitsObject):
    name = "gaussian"
    syntax = {"fwhm": ("tw", Real(), 0),
              "number": ("tn", Integer(), 1),
              "delta": ("td", Real(), 2),
              "cutoff": ("tc", Real(), None)}

    _grammar = r'''start:  "t-type" "=" "2" "\n" normal ~ 3..4
    normal: @assign_among|self.inv_syntax().keys()|
    %ignore SPACE
    '''
    def normal(self, tr):
        return tr[0]

    def tree(self):
        tree = [Tree(Token("RULE", "normal"), [Tree(Token("RULE", "tw"), [Token("computation", self.fwhm)])]),
                Tree(Token("RULE", "normal"), [Tree(Token("RULE", "tn"), [Token("POSINT", self.number)])]),
                Tree(Token("RULE", "normal"), [Tree(Token("RULE", "td"), [Token("computation", self.delta)])])]
        if self.center is not None:
            tree.append(Tree(Token("RULE", "normal"), [Tree(Token("RULE", "tc")), [Token("computation", self.cutoff)]]))

        return Tree(Token("RULE", "start"), tree)

class TimeDistributionBins(PhitsObject):
    name = "timedistbins"
    syntax = {"bins": ("ntt", Integer(), 0),
              "particle_production": ("o-type", Choice10(), None)}

    _grammar = r'''start: t3 | t4
    t3: "t-type" "=" "3" "\n" normal
    t4: "t-type" "=" "4" "\n" normal @assign_then_grid|"particle_production"|
    normal: @assign_then_grid|"bins"|
    %ignore SPACE
    '''

    def o_type(self, tr):
        if len(tr) == 1:
            return ("particle_production", None)
        else:
            return ("particle_production", tr[1][0])

    def ntt(self, tr):
        return ("bins", tr[1])

    def normal(self, tr):
        return tr[0]

    def t3(self, tr):
        return list(tr)

    def t4(self, tr):
        return list(tr)

    def tree(self):
        tree = [Tree(Token("RULE", "normal"), [Tree(Token("RULE", "tw"), [Token("computation", self.width)])]),
                Tree(Token("RULE", "normal"), [Tree(Token("RULE", "tn"), [Token("POSINT", self.number)])]),
                Tree(Token("RULE", "normal"), [Tree(Token("RULE", "td"), [Token("computation", self.delta)])])]
        if self.particle_production is not None:
            tree.append(Tree(Token("RULE", "normal"), [Tree(Token("RULE", "t0")), [Token("computation", self.center)]]))

        return Tree(Token("RULE", "start"), tree)

    




class TimeDistributionFunction(PhitsObject):
    name = "timedistfunction"
    syntax = {"function": ("h(x)", Function(), 0),
              "n_bins": ("ll", Integer(), 1),
              "bounds": (("tg1", "tg2"), (Real(), Real()), 2),
              "particle_production": ("o-type", Choice10(), None)}

    _grammar = r'''start: t5 | t6
    t5: "t-type" "=" "5" "\n" normal ~ 4
    t6: "t-type" "=" "6" "\n" normal ~ 4 @assign_then_grid|"particle_production"|?
    normal: @assign_among|set(self.inv_syntax()) - {"o-type"}|
    %ignore SPACE
    '''

    def o_type(self, tr):
        if len(tr) == 1:
            return ("particle_production", None)
        else:
            return ("particle_production", tr[1][0])

    def assignment(self, tr):
        return (self.inv_syntax()[tr[0]][0], self.inv_syntax()[tr[0]][2].python(tr[1]))

    def normal(self, tr):
        return tr[0]

    def t6(self, tr):
        return list(tr)

    def t5(self, tr):
        return list(tr)

    # TODO: def tree(self):





class TimeDistributionCustom(PhitsObject):
    name = "timedistcustom"
    syntax = {"bounds": (("tg1", "tg2"), (Real(), Real()), 0)}
    _grammar = r'start: "t-type" "=" "100" "\n" @assign_among|["bounds"]|'








class AngleDistributionBins(PhitsObject):
    name = "timedistbins"
    syntax = {"unit": ("a-type", FinBij({"cos": 1, "degree": 11}), 0),
              "bins": ("na", Integer(), 0),
              "particle_production": ("q-type", Choice10(), None)}

    _grammar = r'''start: a1 | a4 | a14 | a11
    a1: "a-type" "=" "1" normal
    a11: "a-type" "=" "11" normal
    a4: "a-type" "=" "4" normal @assign_then_grid|"q-type"|?
    typ: @assign_among|["a-type"]|
    normal: @assign_then_grid|"na"|
    %ignore SPACE
    '''

    def p_type(self, tr):
        if len(tr) == 1:
            return ("particle_production", None)
        else:
            return ("particle_production", tr[1][0])

    def na(self, tr):
        return ("bins", tr[1])

    def normal(self, tr):
        return tr[0]

    def a1(self, tr):
        return list(tr)

    def a11(self, tr):
        return list(tr)






class AngleDistributionFunction(PhitsObject):
    name = "timedistfunction"
    syntax = {"function": ("h(x)", Function(), 0),
              "n_bins": ("ll", Integer(), 1),
              "bounds": (("tg1", "tg2"), (Real(), Real()), 2),
              "particle_production": ("o-type", Choice10(), None)}

    _grammar = r'''start: t5 | t6
    t5: "t-type" "=" "5" "\n" normal ~ 4
    t6: "t-type" "=" "6" "\n" normal ~ 4 @assign_then_grid|self.syntax["particle_production"]|
    normal: @assign_among|{k: self.inv_syntax()[k] for k in set(self.inv_syntax()) - {"o-type"}}|
    %ignore SPACE
    '''

    def o_type(self, tr):
        if len(tr) == 1:
            return ("particle_production", None)
        else:
            return ("particle_production", tr[1][0])

    def assignment(self, tr):
        return (self.inv_syntax()[tr[0]][0], self.inv_syntax()[tr[0]][2].python(tr[1]))

    def normal(self, tr):
        return tr[0]

    def t6(self, tr):
        return list(tr)

    def t5(self, tr):
        return list(tr)

    # TODO: def tree(self):



# class AngleDistribution(PhitsObject):
#     _parser = r"""
#     atype: a1 | a4 | a5 | a6

#     a1:  "a-type"  "="  ("1" | "11")   "\n" assignment numbergrid
#     a4:  "a-type"  "="  ("4" | "14")   "\n" assignment numbergrid assignment numbergrid
#     a5:  "a-type"  "="  ("5" | "15")   "\n" assignment ~ 4
#     a6:  "a-type"  "="  /6|16/   "\n" assignment ~ 2 (assignment numbergrid?)?
#     """

# class EnergyDistribution(PhitsObject):
#     _parser = r"""
#     etype: e1 | e2 | e3 | e4 | e5 | e6 | e7 | e20 | e25 | e28

#     e1:  "e-type"  "="  ("1" | "8" | "11" | "18" | "21" | "22" | "31" | "32")  "\n" assignment numbergrid
#     e4:  "e-type"  "="  ("4" | "9" | "14" | "19" | "23" | "24" | "33" | "34")   "\n" assignment numbergrid assignment numbergrid
#     e2:  "e-type"  "="  ("2" | "12")   "\n" assignment ~ 4
#     e3:  "e-type"  "="  ("3")   "\n" assignment ~
#     e5:  "e-type"  "="  ("5" | "15")   "\n" assignment ~ 4
#     e6:  "e-type"  "="  ("6" | "16")   "\n" assignment ~ 5 numbergrid
#     e7:  "e-type"  "="  ("7")   "\n" assignment ~ 6 numbergrid
#     e20:  "e-type"  "="  ("20")   "\n" assignment
#     e25:  "e-type"  "="  ("25" | "26")   "\n" assignment ~ 14
#     e28:  "e-type"  "="  ("28" | "29")   "\n" assignment ~ 7
#     """
