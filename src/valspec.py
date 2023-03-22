from functools import partial
from hypothesis.strategies import *
# TODO: think about if the python() methods are necessary

class ValSpec():
    def __init__(self, strat):
        self.strat = strat
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
        elif not isinstance(self, OneOf) and not isinstance(other, OneOf):
            return OneOf(self, other)


class Choice10(ValSpec):
    def __init__(self, c_style=False, true=True, false=False):
        super().__init__(one_of(false, true))
        self.c_style = c_style
        self.true = true
        self.false = false

    def phits(self, val):
        if val == self.true:
            return 0 if self.c_style else 1
        elif val == self.false:
            return 1 if self.c_style else 0
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be either True or False; got {va}."), val)

    def python(self, val):
        if val == 0:
                return self.true if self.c_style else self.false
        elif val == 1:
            return self.false if self.c_style else self.true
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be either 0 or 1; got {va}."), val)


class Integer(ValSpec):
    def __init__(self):
        super().__init__(integers())
    def phits(self, val):
        if isinstance(val, int):
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be an integer; got {va}."), val)

    def python(self, val):
        if val % 1 == 0:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be an integer; got {va}."), val)

class Real(ValSpec):
    def __init__(self):
        super().__init__(floats())
    def phits(self, val):
        if isinstance(val, float):
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a float; got {va}."), val)

    def python(self, val):
        return val

class PosInt(ValSpec):
    def __init__(self):
        super().__init__(ints(min_value=0))

    def phits(self, val):
        if isinstance(val, int) and val > 0:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a positive integer; got {va}."), val)

    def python(self, val):
        if val > 0:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be positive; got {va}."), val)




class PosReal(ValSpec):
    def __init__(self):
        super().__init__(floats(min_value=0))
    def phits(self, val):
        if isinstance(val, float) and val > 0:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a positive float; got {va}."), val)

    def python(self, val):
        if val > 0:
            return val
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be positive; got {va}."), val)






class NegDisable(ValSpec):
    def __init__(self):
        super().__init__(one_of(none(), integers(), floats()))
    def phits(self, val):
        if isinstance(val, float) or isinstance(val, int) and val > 0:
            return val
        elif val is None:
            return -1.0
        else:
            return partial(lambda va, var: ValueError(f"`{var}` must be a positive integer or None; got {va}."), val)

    def python(self, val):
        if val > 0:
            return val
        else:
            return None



class Between(ValSpec):
    def __init__(self, start, stop):
        super().__init__(integers(min_vale=start, max_value=stop))
        self.start = start
        self.stop = stop

    def phits(self, val):
        if isinstance(val, int) and val >= self.start and val <= self.stop:
            return val
        else:
            return partial(lambda va, start, stop, var: ValueError(f"`{var}` must be an integer between {start} and {stop}; got {va}."),
                           val, self.start, self.stop)


    def python(self, val):
        if val >= self.start and val <= self.stop:
            return val
        else:
            return partial(lambda va, start, stop, var: ValueError(f"`{var}` must be an integer between {start} and {stop}; got {va}."),
                           val, self.start, self.stop)




class ZeroSpecial(ValSpec):
    def __init__(self, zero):
        super().__init__(one_of(just(zero), integers()))
        self.zero = zero

    def phits(self, val):
        if isinstance(val, int):
            if val == self.zero:
                return 0
            else:
                return val
        else:
            return partial(lambda va, zero, var: ValueError(f"`{var}` must be a positive integer or {zero}; got {va}."), val, self.zero)

    def python(self, val):
        if val == 0:
            return self.zero
        else:
            return val



class FinBij(ValSpec):
    def __init__(self, dic):
        super().__init__(sampled_from(dic.keys()))
        self.dic = dic

    def phits(self, val):
        if val in self.dic:
            return self.dic[val]
        else:
            return partial(lambda va, keys, var: ValueError(f"`{var}` must be one of {keys}; got {va}."), val, list(self.dic.keys()))


    def python(self, val):
        rev = {v: k for k, v in self.dic.items()}
        if val in rev:
            return rev[val]
        else:
            return partial(lambda va, keys, var: ValueError(f"`{var}` must be one of {keys}; got {va}."), val, list(rev.keys()))

class IsA(ValSpec):
    def __init__(self, cls, index=False):
        super().__init__(from_type(cls)) # TODO: sus
        self.cls = cls
        self.index = index

    def phits(self, val):
        if not isinstance(val, self.cls):
            return partial(lambda va, cls, var: ValueError(f"`{var}` must be an instance of {cls}; got {val}."), val, self.cls)

        if self.index:
            return val.index
        else:
            return val.definition()

    def python(self, val):
        return val


# TODO: needed?
class OneOf(ValSpec):
    def __init__(self, *args):
        assert all(map(lambda x: isinstance(x, ValSpec), args)), "All arguments to OneOf must be value specifications."
        super().__init__(one_of(*[i.strat for i in args]))
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
