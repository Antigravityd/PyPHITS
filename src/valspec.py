
class ValSpec():
    def __init__(self, parse_rule):
        self.parse_rule = parse_rule
    def __or__(self, other):
        # TODO: think about if copying is necessary
        if isinstance(self, OneOf) and isinstance(other, OneOf):
            self.choices.append(other.choices)
            return self
        elif isinstance(self, OneOf) and not isinstance(other, OneOf):
            self.choices.append(other)
            return self
        elif not isinstance(self, OneOf) and isinstance(other, OneOf):
            other.choices.append(self)
            return other
        elif not isinstance(self, OneOf) and not isinstance(other, OneOf)
        return OneOf(self, other)


class Choice10(ValSpec):
    def __init__(self, c_style=False, true=True, false=False):
        super().__init__("INT")
        self.c_style = c_style
        self.true = true
        self.false = false

    def phits(self, val):
        if val == self.true:
            return 0 if self.c_style else 1
        elif val == false:
            return 1 if c_style else 0
        else:
            return lambda var: ValueError(f"`{var}` must be either True or False; got {val}.")

    def python(self, val):
        if val == 0:
                return true if c_style else false
        elif val == 1:
            return false if c_style else true
        else:
            return lambda var: ValueError(f"`{self.ident_map[var]}` must be either 0 or 1; got {val}.")


class Integer(ValSpec):
    def __init__(self):
        super().__init__("INT")

    def phits(self, val):
        if isinstance(val, int):
            return val
        else:
            return lambda var: ValueError(f"`{var}` must be an integer; got {val}.")

    def python(self, val):
        if x % 1 == 0:
            return val
        else:
            return lambda var: ValueError(f"`{self.ident_map[var]}` must be an integer; got {val}.")

class Real(ValSpec):
    def __init__(self):
        super().__init__("NUMBER")

    def phits(self, val):
        if isinstance(val, float):
            return val
        else:
            return lambda var: ValueError(f"`{var}` must be a float; got {val}.")

    def python(self, val):
        return val

class PosInt(ValSpec):
    def __init__(self):
        super().__init__("POSINT")

    def phits(self, val):
        if isinstance(val, int) and val > 0:
            return val
        else:
            return lambda var: ValueError(f"`{var}` must be a positive integer; got {val}.")

    def python(self, val):
        if val > 0:
            return val
        else:
            return lambda var: ValueError(f"`{self.ident_map[var]}` must be positive; got {val}.")




class PosReal(ValSpec):
    def __init__(self):
        super().__init__("NUMBER")
    def phits(self, val):
        if isinstance(val, float) and val > 0:
            return val
        else:
            lambda var: ValueError(f"`{var}` must be a positive floating-point value; got {val}.")

    def python(self, val):
        if val > 0:
            return val
        else:
            lambda var: ValueError(f"`{self.ident_map[var]}` must be positive; got {val}.")




# incorrect
# class TripleChoice():
#     parse_rule = "INT"

#     def produce(self):
#         return tokmap("INT")

#     def phits(self):
#         if isinstance(val, int):
#             return True if val > 0 else False
#         else:
#             return lambda var: ValueError(f"`{var}` must be an integer; got {val}.")

#     def python(self):
#         if val > 0:
#             return True
#         else:
#             return False







class NegDisable(ValSpec):
    def __init__(self):
        super().__init__("NUMBER")
    def phits(self, val):
        if isinstance(val, float) or isinstance(val, int) and val > 0:
            return val
        elif val is None:
            return -1.0
        else:
            # TODO: think about how to make this work
            return lambda s, var: ValueError(f"`{s.ident_map[var]}` must be a positive integer or None; got {val}.")

    def python(self, val):
        if val > 0:
            return val
        else:
            return None



class Between(ValSpec):
    def __init__(self, start, stop):
        super().__init__("INT")
        self.start = start
        self.stop = stop

    def phits(self, val):
        if isinstance(val, int) and val >= self.start and val <= self.stop:
            return val
        else:
            return lambda var: ValueError(f"`{var}` must be an integer between {self.start} and {self.stop}, inclusive; got {val}.")

    def python(self, val):
        if val >= self.start and val <= self.stop:
            return val
        else:
            return lambda var: ValueError(f"`{self.ident_map[var]}` must be between {self.start} and {self.stop}, inclusive; got {val}.")



class ZeroSpecial(ValSpec):
    def __init__(self, zero):
        super().__init__("INT")
        self.zero = zero

    def phits(self, val):
        if isinstance(val, int):
            if val == 0:
                return self.zero
            else:
                return val
        else:
            return lambda var: ValueError(f"`{var}` must be an integer; got {val}.")

    def python(self, val):
        if val == zero:
            return 0
        else:
            return val



class FinBij(ValSpec):
    def __init__(self, dic):
        super().__init__("INT") # TODO: is this all that's needed?
        self.dic = dic

    def phits(self, val):
        if val in self.dic:
            return self.dic[val]
        else:
            return lambda var: ValueError(f"`{var}` must be one of {list(dic.keys())}; got {val}.")

    def python(self, val):
        rev = {v: k for k, v in self.dic.items()}
        if val in rev:
            return rev[val]
        else:
            return lambda var: ValueError(f"`{self.ident_map[var]}` must be one of {list(rev.keys())}; got {val}.")

class IsA(ValSpec):
    def __init__(self, typ, index=False):
        super.__init__(f"@parseof|{typ}|") # TODO: sus
        self.index = index
        assert type(typ) == type, "Pass class name, not instance, to IsA."
        assert issubclass(typ, PhitsObject), "IsA must be passed a PhitsObject."
        self.typ = typ

    def phits(self, val):
        if self.index:
            return val.index
        else:
            return val.definition()

    def python(self, val):
        if self.index:
            if isinstance(val, int):
                return val
            else:
                return lambda var: ValueError(f"`{var}` must be an integer; got {val}.")
        else:
            if isinstance(val, typ):
                return val
            else:
                return lambda var: ValueError(f"`{var}` must be an object of type {typ}; got {val}.")


class OneOf(ValSpec):
    def __init__(self, *args):
        assert all(map(lambda x: isinstance(x, ValSpec), args)), "All arguments to OneOf must be value specifications."
        super().__init__("(" + " | ".join([i.parse_rule for i in args]) + ")")
        self.choices = args

    def phits(self, val):
        return self.that_which_applies(val).phits(val)

    def python(self, val):
        return self.that_which_applies(val).python(val)

    def that_which_applies(self, val):
        def _applies(s, val):
            try:
                s.phits(val)
                return ("phits", s)
            except Exception:
                try:
                    s.python(val)
                    return ("python", s)
                except Exception:
                    return None

        applicable = filter(lambda x: _applies(x, val) is not None, self.choices)
        if len(applicable) == 0:
            return ValueError("Empty OneOf value specification.")
        else:
            assert len(applicable) == 1, "Ambiguous OneOf value specification."
            return applicable[0]
