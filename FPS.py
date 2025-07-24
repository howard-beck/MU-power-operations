import math
from sympy import *
import sympy
from utils import *

alpha = symbols("alpha")
beta = symbols("beta")

class MultiIndex:
    def __init__(self, idx):
        if type(idx) is MultiIndex:
            idx = idx.idx
        assert type(idx) is dict

        self.idx = idx

        self.dims = len(self.idx)
        
    def combine_vars(self, other):
        if type(other) is MultiIndex:
            return self.combine_vars(self, other.idx)
        
        return set(self.idx) | set(other)
    
    def __lt__(self, midx):
        if type(midx) is MultiIndex:
            return self < midx.idx
        assert type(midx) is dict

        vars = self.combine_vars(midx)

        for var in vars:
            if var in self.idx and var not in midx:
                # can't be less, this coordinate is more
                if self.idx[var] > 0:
                    return False
            if var in self.idx and var in midx:
                if self.idx[var] < midx[var]:
                    return True
        return False
    
    def __eq__(self, midx):
        if type(midx) is MultiIndex:
            return self < midx.idx
        assert type(midx) is dict

        vars = self.combine_vars(midx)

        for var in vars:
            if var in self.idx and var in midx:
                if self.idx[var] != midx[var]:
                    return False
            if var not in self.idx and var in midx:
                if midx[var] != 0:
                    return False
            if var in self.idx and var not in midx:
                if self.idx[var] != 0:
                    return False
        return False
    
    def __gt__(self, midx):
        if type(midx) is MultiIndex:
            return self > midx.idx
        assert type(midx) is dict

        vars = self.combine_vars(midx)

        for var in vars:
            if var in midx and var not in self.idx:
                # can't be more, this coordinate is less
                if midx[var] > 0:
                    return False
            if var in self.idx and var in midx:
                if self.idx[var] < midx[var]:
                    return True
        return False
    
    def __le__(self, midx):
        return (self < midx) or (self == midx)
    def __ge__(self, midx):
        return (self > midx) or (self == midx)

    def __add__(self, midx):
        if type(midx) is MultiIndex:
            return self + midx.idx
        
        combine_dict = {var: self.idx[var] for var in self.idx}

        for var in midx:
            if var in combine_dict:
                combine_dict[var] += midx[var]
            else:
                combine_dict[var] = midx[var]
        
        return MultiIndex(combine_dict)

    def max(self, midx):
        if type(midx) is MultiIndex:
            return self.max(midx.idx)
        
        combine_dict = {var: self.idx[var] for var in self.idx}

        for var in midx:
            if var in self.idx:
                combine_dict[var] = max(midx[var], self.idx[var])
            else:
                combine_dict[var] = midx[var]
        
        return MultiIndex(combine_dict)


class FPS:
    def __init__(self, generator, vars, name = None, max_degree = None):
        self.generator = generator
        self.name = name
        
        if type(vars) is Symbol:
            vars = {vars}
        self.vars = set(vars)
        self.var_names = [var.name for var in vars]
        self.max_degree = max_degree
        if max_degree is None:
            self.max_degree = MultiIndex({
                var.name: float("inf")
                for var in vars
            })
        self.dims = len(self.vars)

        self.coeffs = {}

        # sometimes, can be useful to save the powers
        self.powers = {}
        self.save_powers = False

        self.is_symmetric = False

        self.on_access = None
        self.on_calculate = None
    
    def as_parseable_obj(self):
        obj = {
            "val": {}
        }
        for idx in self.coeffs:
            obj["val"][",".join(str(i) for i in idx)] = str(self.coeffs[idx])
        if len(self.powers) > 0:
            obj["powers"] = {}
            for power in self.powers:
                obj["powers"][str(power)] = self.powers[power].as_parseable_obj()
        obj["vars"] = self.var_names
        return obj
    
    def from_parseable_obj(self, obj):
        assert "val" in obj
        for idx in obj["val"]:
            new_idx = tuple(int(i) for i in idx.split(","))
            self.coeffs[new_idx] = parse_expr(
                obj["val"][idx],
                {
                    "alpha": alpha,
                    "beta": beta
                }
            )
        if "powers" in obj:
            self.save_powers = True
            for power in obj["powers"]:
                if power not in self.powers:
                    self.powers[int(power)] = self ** int(power)
                self.powers[int(power)].from_parseable_obj(obj["powers"][power])
        if "vars" in obj:
            self.var_names = obj["vars"]

    
    def coeff(self, idx):
        idx_tuple = None
        if type(idx) is int:
            assert self.dims == 1
            idx_tuple = (idx,)
        elif type(idx) is tuple:
            idx_tuple = idx
            assert len(idx) == self.dims
        else:
            assert type(idx) is dict

            idx_dict = {}
            for var_name in idx:
                if var_name not in self.var_names:
                    if idx[var_name] != 0:
                        return 0
            for var_name in self.var_names:
                if var_name not in idx:
                    idx_dict[var_name] = 0
                else:
                    idx_dict[var_name] = idx[var_name]
            idx_tuple = tuple(idx_dict[var_name] for var_name in self.var_names)

        if self.on_access is not None:
            self.on_access(idx_tuple)
        
        if idx_tuple in self.coeffs:
            return self.coeffs[idx_tuple]
        if self.is_symmetric:
            assert self.dims == 2
            if idx_tuple[1] > idx_tuple[0]:
                return self.coeff((idx_tuple[1], idx_tuple[0]))
        
        if self.dims == 1:
            try:
                self.coeffs[idx_tuple] = self.generator(idx_tuple)
            except TypeError:
                self.coeffs[idx_tuple] = self.generator(idx_tuple[0])
        else:
            self.coeffs[idx_tuple] = self.generator(idx_tuple)
        if self.on_calculate is not None:
            self.on_calculate(idx_tuple)
        return self.coeffs[idx_tuple]
    
    def __add__(self, g):
        new_vars = list(self.vars | g.vars)

        def generator(midx):
            idx = {
                new_vars[i].name: midx[i]
                for i in range(len(new_vars))
            }

            if self.max_degree < idx:
                return g.coeff(idx)
            elif g.max_degree < idx:
                return self.coeff(idx)
            
            return self.coeff(idx) + g.coeff(idx)
    
        return FPS(
            lambda idx: generator(idx),
            vars = self.vars | g.vars,
            max_degree = self.max_degree.max(g.max_degree)
        )
    
    def __sub__(self, g):
        new_vars = list(self.vars | g.vars)

        def generator(midx):
            idx = {
                new_vars[i].name: midx[i]
                for i in range(len(new_vars))
            }

            if self.max_degree < idx:
                return -g.coeff(idx)
            elif g.max_degree < idx:
                return self.coeff(idx)
            
            return self.coeff(idx) - g.coeff(idx)
    
        return FPS(
            lambda idx: generator(idx),
            vars = self.vars | g.vars,
            max_degree = self.max_degree.max(g.max_degree)
        )
    
    def __mul__(self, g):
        if type(g) is int or \
            type(g) is float or \
            type(g) is sympy.core.numbers.Integer or \
            type(g) is sympy.core.numbers.Float or \
            type(g) is sympy.core.numbers.Rational or \
            type(g) is sympy.core.symbol.Symbol or \
            type(g) is sympy.core.numbers.One or \
            type(g) is sympy.core.power.Pow or \
            type(g) is sympy.core.numbers.Half:
                # don't use map, since this lets you easily compute powers:
                return ScalarMultiple(self, g)

        total_vars = list(self.vars | g.vars)

        def generator(idx):
            indices_to_eval = []

            def index_recursion(partial_idx):
                m = len(partial_idx)
                if m == len(total_vars):
                    indices_to_eval.append({
                        p[0]: p[1]
                        for p in partial_idx
                    })
                    return
                for i in range(idx[m]+1):
                    var_name = total_vars[m].name
                    index_recursion(partial_idx | {(var_name, i)})
            index_recursion(set())

            ret = 0
            for midx in indices_to_eval:
                oidx = {
                    total_vars[i].name: idx[i] - midx[total_vars[i].name]
                    for i in range(len(total_vars))
                }
                a = self.coeff(midx)
                b = g.coeff(oidx)

                ret += a * b
            return ret

        return FPS(
            generator,
            vars = self.vars | g.vars,
            max_degree = self.max_degree + g.max_degree
        )

    def __pow__(self, n):
        assert n >= 0
        assert type(n) is int


        if n == 0:
            return FPS.get_const(self.vars, 1)
        if n == 1:
            return self
        if self.save_powers:
            if n in self.powers:
                return self.powers[n]
        
        a = n//2
        b = n - a
        ret = (self**a) * (self**b)

        if self.save_powers:
            self.powers[n] = ret
        return ret
    
    def comp(self, assignments, save_terms = False):
        return CompositeFPS(self, assignments, save_terms)

    def comp_inv(self):
        assert self.dims == 1

        assert self.coeff(0) == 0 or self.coeff(0) == 0.0
        assert self.coeff(1) == 1 or self.coeff(1) == 1.0

        def generator(n):
            if n == 0:
                return 0
            if n == 1:
                return 1
            
            derivs = [self.coeff(i) * factorial(i) for i in range(n+1)]
            f_hat = [Rational(1, i) * derivs[i] for i in range(2, n+1)]

            ret = 0

            for k in range(1, n):
                rising = 1
                for i in range(0, k):
                    rising *= (n+i)
                ret += (-1)**k * rising * bell(n-1, k, f_hat)

            return simplify(Rational(1, factorial(n)) * ret)
        return FPS(generator, self.vars)

    def mult_inv(self, var = None, const_inv = None):
        assert self.dims == 1

        return MultiplicativeInverse1D(
            self,
            self.vars,
            const_inv
        )
    
    def shift(self, amount):
        assert self.dims == 1

        def generator(idx):
            return self.coeff(idx - amount)
        return FPS(
            generator,
            self.vars
        )
    
    def map_coeffs(self, f):
        def generator(idx):
            return f(self.coeff(idx))
        return FPS(generator, self.vars)



    def print(self, order = 2):
        indices_to_eval = []
        def index_recursion(partial_idx):
            total_partial_index = sum(partial_idx)
            if len(partial_idx) == self.dims:
                indices_to_eval.append(tuple(partial_idx))
            else:
                for i in range(order - total_partial_index + 1):
                    index_recursion(partial_idx + [i])
        index_recursion([])

        if self.name is not None:
            print(self.name + " = ")
        for idx in indices_to_eval:
            print(str(idx) + ": " + str(self.coeff(idx)))
    
    def calculate_up_to(self, order = 5, callback = lambda: (), every = 10):
        def get_multi_indices_total_n(n):
            indices = []
            def recursor(partial_idx):
                total = sum(partial_idx)
                if len(partial_idx) == self.dims - 1:
                    indices.append(tuple([n - total] + partial_idx))
                else:
                    for i in range(n - total+1):
                        recursor([i] + partial_idx)
            recursor([])
            return indices

        orig_on_calculate = self.on_calculate

        i = 0
        def new_on_calculate(idx):
            nonlocal i

            i += 1
            if i % every == 0:
                callback()
            if orig_on_calculate is not None:
                orig_on_calculate(idx)
        self.on_calculate = new_on_calculate

        if self.name is not None:
            print("Calculating " + self.name)
        for n in range(order+1):
            for idx in get_multi_indices_total_n(n):
                print(str(idx) + ": " + str(self.coeff(idx)))
        
        self.on_calculate = orig_on_calculate
        callback()
    
    """
    def mult_fps_generator(f, g, idx):
        if f.max_degree + g.max_degree < idx:
            return 0
        
        ret = 0
        for i in range()
        return sum(
            f.coeff(i) * g.coeff(idx - i)
            for i in range(idx+1)
        )
    """

    def get_one_term(vars, coeff, idx):
        vars = list(vars)
        idx_tuple = [(vars[i].name, idx[i]) for i in range(len(vars))]
        def generator(oidx):
            if oidx == idx:
                return coeff
            return 0
        return FPS(
            generator,
            vars = vars,
            max_degree = MultiIndex({
                vars[i].name: idx[i]
                for i in range(len(vars))
            })
        )
    
    def get_const(vars, c):
        if type(vars) is Symbol:
            return FPS.get_one_term([vars], c, (0,))
        else:
            return FPS.get_one_term(vars, c, tuple([0]*len(vars)))
    

    def from_polynomial(p, vars):
        if type(vars) is Symbol:
            vars = [vars]

        def generator(idx):
            ret = p
            for i in range(len(vars)):
                ret = ret.coeff(vars[i], idx[i])
            return ret

        return FPS(
            generator,
            vars
        )
    
    def prod(fs):
        prod = None

        for f in fs:
            if prod is None:
                prod = f
            else:
                prod *= f
        assert prod is not None
        return prod

class CompositeFPS(FPS):
    def __init__(self, f, assignments, save_terms = False):
        self.f = f

        assert len(assignments) == self.f.dims

        vars_to_change = [a[0].name for a in assignments]
        vals_to_change = [a[1] for a in assignments]

        new_vars = set()
        for val in vals_to_change:
            new_vars |= val.vars

            tuple0 = tuple([0]*len(val.var_names))
            coeff0 = val.coeff(tuple0)

            assert coeff0 == 0 or coeff0 == 0.0

        assert set(vars_to_change) == set(self.f.var_names)

        self.map = {
            vars_to_change[i]: vals_to_change[i]
            for i in range(len(assignments))
        }

        self.cross_powers = {}
        self.save_terms = save_terms

        super().__init__(self.generator, new_vars)
    
    def get_cross_terms(self, midx):
        assert all([i >= 0 for i in midx])

        if self.save_terms and midx in self.cross_powers:
            return self.cross_powers[midx]
        ret = None
        if not self.save_terms or midx not in self.cross_powers:
            prod = FPS.get_const(self.vars, 1)
            for i in range(len(self.f.var_names)):
                prod *= self.map[self.f.var_names[i]]**midx[i]
            ret = prod
        if self.save_terms:
            self.cross_powers[midx] = ret
        return ret



    def generator(self, idx):
        total = 0
        if self.f.dims == 1 and type(idx) is int:
            total = idx
        else:
            total = sum(idx)

        indices_to_eval = set()
        def index_recursion(partial_idx):
            total_partial_index = sum(partial_idx)
            if len(partial_idx) == len(self.f.var_names):
                indices_to_eval.add(tuple(partial_idx))
            else:
                for i in range(total - total_partial_index + 1):
                    index_recursion(partial_idx + [i])
        index_recursion([])

        ret = 0
        for midx in indices_to_eval:
            a = self.f.coeff(midx)
            prod = self.get_cross_terms(midx)
            ret += a * prod.coeff(idx)

        return simplify(ret)

class ScalarMultiple(FPS):
    def __init__(self, f, k):
        super().__init__(self.generator, f.vars)

        self.f = f
        self.k = k
    
    def generator(self, idx):
        return self.k * self.f.coeff(idx)
    
    def __pow__(self, n):
        return (self.f**n) * (self.k**n)

class IteratedMultivariableFPS(FPS):
    def __init__(self, term_generator, vars, out_var):
        self.terms = {}
        
        self.term_generator = term_generator
        self.out_var = out_var

        super().__init__(self.term_generator, vars)

    def term(self, i):
        if i not in self.terms:
            the_term = self.term_generator(i)
            assert the_term.vars == self.vars - {self.out_var}
            self.terms[i] = the_term
        return self.terms[i]
    
    def term_generator(self, idx):
        in_dict = {
            self.vars[i]: idx[i]
            for i in range(len(idx))
            if self.vars[i] != self.out_var
        }

        out_idx = idx[self.var_names.index(self.out_var.name)]

        the_term = self.term(out_idx)
        return the_term.coeff(in_dict)


class MultiplicativeInverse1D(FPS):
    def __init__(self, f, vars, const_inv = None):
        assert f.dims == 1

        if const_inv is None:
            const_inv = f.coeff(0)**(-1)
        self.const_inv = const_inv
        self.f = f

        super().__init__(self.term_generator, vars)

    def term_generator(self, n):
        if n == 0:
            return self.const_inv
        else:
            ret = 0
            for i in range(n):
                ret += self.coeff(i) * self.f.coeff(n - i)
            return ret * self.const_inv

class MultiplicativeInverseND(FPS):
    def __init__(self, f, vars, out_var, const_inv):
        super().__init__(self.term_generator, vars)

    def term_generator(self, n):
        if n == 0:
            return self.const_inv
        else:
            ret = 0
            for i in range(n):
                ret += self.coeff(i) * self.f.coeff(n - i)
            return ret * self.const_inv