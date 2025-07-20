import math
from sympy import *
from utils import *

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
    
    def get_multi_index(self, idx):
        if self.dims == 1 and type(idx) is int:
            return {
                self.var_names[0]: idx
            }
        
        return {
            self.var_names[i]: idx[i]
            for i in range(len(self.var_names))
        }
    
    def get_coeff(self, idx):
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

        if idx_tuple in self.coeffs:
            return self.coeffs[idx_tuple]
        
        if self.dims == 1:
            try:
                self.coeffs[idx_tuple] = self.generator(idx_tuple)
            except TypeError:
                self.coeffs[idx_tuple] = self.generator(idx_tuple[0])
        else:
            self.coeffs[idx_tuple] = self.generator(idx_tuple)
        return self.coeffs[idx_tuple]
    
    def __add__(self, g):
        return FPS(
            lambda idx: FPS.add_fps_generator(self, g, self.get_multi_index(idx)),
            vars = self.vars | g.vars,
            max_degree = self.max_degree.max(g.max_degree)
        )
    
    def __sub__(self, g):
        return FPS(
            lambda idx: FPS.sub_fps_generator(self, g, self.get_multi_index(idx)),
            vars = self.vars | g.vars,
            max_degree = self.max_degree.max(g.max_degree)
        )
    
    def __mul__(self, g):
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
                a = self.get_coeff(midx)
                b = g.get_coeff(oidx)

                ret += a * b
            return ret

        return FPS(
            generator,
            vars = self.vars | g.vars,
            max_degree = self.max_degree + g.max_degree
        )

    def __pow__(self, n):
        assert n >= 0

        if n == 0:
            return FPS.get_const(self.vars, 1)
        if n == 1:
            return self
        else:
            return self**(n - 1) * self
    
    def comp(self, assignments):
        assert len(assignments) == self.dims

        vars_to_change = [a[0].name for a in assignments]
        vals_to_change = [a[1] for a in assignments]

        assert set(vars_to_change) == set(self.var_names)

        map = {
            vars_to_change[i]: vals_to_change[i]
            for i in range(len(assignments))
        }

        new_vars = set()
        for val in vals_to_change:
            new_vars |= val.vars

            tuple0 = tuple([0]*len(val.var_names))
            coeff0 = val.get_coeff(tuple0)

            assert coeff0 == 0 or coeff0 == 0.0
        
        def generator(idx):
            if self.dims == 1 and type(idx) is int:
                total = idx
            else:
                total = sum(idx)

            indices_to_eval = set()
            def index_recursion(partial_idx):
                total_partial_index = sum(partial_idx)
                if len(partial_idx) == len(self.var_names):
                    indices_to_eval.add(tuple(partial_idx))
                else:
                    for i in range(total - total_partial_index + 1):
                        index_recursion(partial_idx + [i])
            index_recursion([])

            ret = 0
            for midx in indices_to_eval:
                a = self.get_coeff(midx)
                prod = FPS.get_const(new_vars, 1)
                for i in range(len(self.var_names)):
                    prod *= (map[self.var_names[i]] ** midx[i])
                ret += a * prod.get_coeff(idx)

            return simplify(ret)

        return FPS(
            generator,
            new_vars
        )



    def comp_inv(self):
        assert self.dims == 1

        assert self.get_coeff(0) == 0 or self.get_coeff(0) == 0.0
        assert self.get_coeff(1) == 1 or self.get_coeff(1) == 1.0

        def generator(n):
            if n == 0:
                return 0
            if n == 1:
                return 1
            
            derivs = [self.get_coeff(i) * factorial(i) for i in range(n+1)]
            f_hat = [Rational(1, i) * derivs[i] for i in range(2, n+1)]

            ret = 0

            for k in range(1, n):
                rising = 1
                for i in range(0, k):
                    rising *= (n+i)
                ret += (-1)**k * rising * bell(n-1, k, f_hat)

            return simplify(Rational(1, factorial(n)) * ret)
        return FPS(generator, self.vars)

    def print(self, order = 2):
        indices_to_eval = set()
        def index_recursion(partial_idx):
            total_partial_index = sum(partial_idx)
            if len(partial_idx) == self.dims:
                indices_to_eval.add(tuple(partial_idx))
            else:
                for i in range(order - total_partial_index + 1):
                    index_recursion(partial_idx + [i])
        index_recursion([])

        for idx in indices_to_eval:
            print(str(idx) + ": " + str(self.get_coeff(idx)))

    def add_fps_generator(f, g, idx):
        if f.max_degree < idx:
            return g.get_coeff(idx)
        elif g.max_degree < idx:
            return f.get_coeff(idx)
        
        return f.get_coeff(idx) + g.get_coeff(idx)
    
    def sub_fps_generator(f, g, idx):
        if f.max_degree < idx:
            return -g.get_coeff(idx)
        elif g.max_degree < idx:
            return f.get_coeff(idx)
        
        return f.get_coeff(idx) - g.get_coeff(idx)
    
    """
    def mult_fps_generator(f, g, idx):
        if f.max_degree + g.max_degree < idx:
            return 0
        
        ret = 0
        for i in range()
        return sum(
            f.get_coeff(i) * g.get_coeff(idx - i)
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
            return FPS.get_one_term(vars, c, (0,))
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
            