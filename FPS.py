import math
from sympy import *
from utils import *

class MultiIndex:
    def __init__(self, idx):
        if type(idx) is MultiIndex:
            idx = idx.idx
        
        self.idx = idx

        self.dims = 1
        if type(idx) is Symbol:
            self.multivariate = False
        else:
            self.multivariate = True
            self.dims = len(idx)
            
        self.ordering = "min"
        

    def __lt__(self, midx):
        if type(midx) is MultiIndex:
            return self < midx.idx
        elif type(midx) is int:
            assert self.multivariate is False
            return self.idx < midx
        assert type(midx) is tuple or type(midx) is list

        assert self.ordering == "lex" or self.ordering == "total" or self.ordering == "min"

        if self.ordering == "lex":
            for i in range(self.dims):
                if self.idx[i] < midx[i]:
                    return True
                elif self.idx[i] > midx[i]:
                    return False
            return False
        elif self.ordering == "total":
            return sum(self.idx) < sum(midx)
        elif self.ordering == "min":
            for i in range(self.dims):
                if self.idx[i] < midx[i]:
                    return True
            return False
    
    def __eq__(self, midx):
        if type(midx) is MultiIndex:
            return self == midx.idx
        elif type(midx) is int:
            assert self.multivariate is False
            return self.idx == midx
        assert type(midx) is tuple or type(midx) is list

        assert self.ordering == "lex" or self.ordering == "total" or self.ordering == "min"

        if self.ordering == "lex":
            return self.idx == midx
        if self.ordering == "total":
            # NOTE: only a partial ordering
            return self.idx == midx
        elif self.ordering == "min":
            return self.idx == midx
    
    def __gt__(self, midx):
        if type(midx) is MultiIndex:
            return self > midx.idx
        elif type(midx) is int:
            assert self.multivariate is False
            return self.idx > midx
        assert type(midx) is tuple or type(midx) is list

        assert self.ordering == "lex" or self.ordering == "total" or self.ordering == "min"

        if self.ordering == "lex":
            for i in range(self.dims):
                if self.idx[i] > midx[i]:
                    return True
                elif self.idx[i] < midx[i]:
                    return False
            return False
        elif self.ordering == "total":
            return sum(self.idx) > sum(midx)
        elif self.ordering == "min":
            for i in range(self.dims):
                if self.idx[i] > midx[i]:
                    return True
            return False
    
    def __le__(self, midx):
        return (self < midx) or (self == midx)
    def __ge__(self, midx):
        return (self > midx) or (self == midx)

    def __add__(self, midx):
        if type(midx) is MultiIndex:
            return self + midx.idx
        if type(midx) is int:
            assert self.multivariate is False
            return MultiIndex(self.idx + midx)
        
        assert self.dims == len(midx)

        return MultiIndex([
            self.idx[i] + midx[i]
            for i in range(len(midx))
        ])
        


class FPS:
    def __init__(self, generator, vars, name = None, max_degree = None):
        self.generator = generator
        self.name = name
        
        self.var = None
        self.vars = []
        self.multivariate = False
        self.max_degree = max_degree
        self.dims = 1
        if type(vars) is Symbol:
            self.var = vars
            self.vars = [vars]
            if max_degree is not None:
                self.max_degree = MultiIndex(float("inf"))
        else:
            self.vars = vars
            self.multivariate = True
            self.dims = len(self.vars)
            if max_degree is not None:
                self.max_degree = MultiIndex([float("inf")]*self.dims)

        self.coeffs = {}

        
    
    def get_coeff(self, idx):
        if idx in self.coeffs:
            return self.coeffs[idx]
        
        self.coeffs[idx] = self.generator(idx)
        return self.coeffs[idx]
    
    def __add__(self, g):
        assert self.vars == g.vars

        return FPS(
            lambda idx: FPS.add_fps_generator(self, g, idx),
            vars = self.vars,
            max_degree = max(self.max_degree, g.max_degree)
        )
    
    def __sub__(self, g):
        assert self.vars == g.vars

        return FPS(
            lambda idx: FPS.sub_fps_generator(self, g, idx),
            vars = self.vars,
            max_degree = max(self.max_degree, g.max_degree)
        )
    
    def __mult__(self, g):
        assert self.vars == g.vars

        return FPS(
            lambda idx: FPS.mult_fps_generator(self, g, idx),
            vars = self.vars,
            max_degree = self.max_degree + g.max_degree
        )

    def __pow__(self, n):
        assert n >= 0

        if n == 0:
            return FPS.get_const(self.vars, 1)
        if n == 1:
            return self
        else:
            return self.pow(n - 1) * self

    def comp_inv(self):
        assert self.multivariate is False

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
    
    def mult_fps_generator(f, g, idx):
        assert f.dims == 1
        assert g.dims == 1

        if f.max_degree + g.max_degree < idx:
            return 0
        
        return sum(
            f.get_coeff(i) * g.get_coeff(idx - i)
            for i in range(idx+1)
        )
    
    def get_one_term(vars, coeff, idx):
        def generator(oidx):
            if oidx == idx:
                return coeff
            return 0
        return FPS(
            generator,
            vars = vars,
            max_degree = MultiIndex(idx)
        )
    
    def get_const(vars, c):
        if type(vars) is Symbol:
            return FPS.get_one_term(vars, c, 0)
        else:
            return FPS.get_one_term(vars, c, [0]*len(vars))