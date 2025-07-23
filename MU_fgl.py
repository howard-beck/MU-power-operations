from sympy import *
from FPS import *
from fgls_auto import *
import time


primes = [2, 3, 5, 7, 11, 13, 17]
last_primality_checked = primes[-1]
def check_primes_until(n):
    global last_primality_checked
    if last_primality_checked >= n:
        return
    
    for j in range(last_primality_checked + 1, n+1):
        is_prime = True
        for p in primes:
            if j % p == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(j)
    new_primes = [p for p in primes if p > last_primality_checked]
    if len(new_primes) > 0:
        print("Added primes: " + str(new_primes))
    last_primality_checked = n

alpha = symbols("alpha")



class ComplexCobordism(FGL2):
    def __init__(self, use_hazewinkel = True):
        self.xs = {}
        self.vs = {}
        self.cs = {}

        self.use_hazewinkel = use_hazewinkel

        super().__init__(
            log_fgl = FPS(
                self.log_fgl_generator,
                vars = alpha,
                name = "log"
            ),
            name = "MU2"
        )

        self.chis = {}

        self.a_series = {}
        self.a = {}

    def v(n):
        if n == 1:
            return 1

        check_primes_until(n)
        
        for p in primes:
            if n % p == 0:
                power = math.log(n) / math.log(p)
                if abs(power - round(power)) <= 0.000001:
                    return p
                else:
                    return 1
        return 1

    def c(p, d):
        vd = ComplexCobordism.v(d)
        if vd == 1 or vd == p:
            return 1
        ret = vd
        while ret % p != 1:
            ret += vd
        
        # = 0 mod vd
        # = 1 mod p
        return ret

    def mu(n, d):
        ret = 1

        check_primes_until(n)

        for p in primes:
            if n % p == 0:
                ret *= ComplexCobordism.c(p, d)
        return ret



    def get_x(self, i):
        if i == 0:
            return 1
        if i not in self.xs:
            self.xs[i] = symbols("x" + str(i))
        return self.xs[i]

    def get_v(self, i):
        if i not in self.vs:
            self.vs[i] = symbols("v" + str(i))
        return self.vs[i]

    def get_c(self, i):
        if i not in self.cs:
            if i == 0:
                self.cs[0] = 1
            else:
                if self.use_hazewinkel:
                    k = i + 1

                    ret = Rational(k, ComplexCobordism.v(k)) * self.get_x(k-1)

                    for d in range(2, k):
                        if k % d == 0:
                            coeff = Rational(ComplexCobordism.mu(k, d) * d, ComplexCobordism.v(d))
                            ret += coeff * self.get_x(d-1)**(round(k/d)) * self.get_c(round(k/d)-1)
                    
                    self.cs[i] = ret
                else:
                    self.cs[i] = symbols("c" + str(i))
        return self.cs[i]

    def log_fgl_generator(self, n):
        if n == 0:
            return 0
        if n == 1:
            return 1
        return self.get_c(n-1)/n

    def chi(self, p):
        if p in self.chis:
            return self.chis[p]
        else:
            self.chis[p] = FPS.prod(
                MU.get_n_series(i)
                for i in range(1, p)
            )
            self.chis[p].name = "chi_" + str(p)
            self.chis[p].save_powers = True
            return self.chis[p]
    
    def get_a_series(self, p):
        if p in self.a_series:
            return self.a_series[p]
        ret = FPS.prod(MU.get_beta_plus_n_series(i) for i in range(1, p))
        ret.name = "a-series @ " + str(p)
        self.a_series[p] = ret
        return ret
    
    def get_a(self, p, i):
        if (p, i) in self.a:
            return self.a[(p, i)]

        a_series = self.get_a_series(p)
        if False:
            self.a[(p, 0)] = self.chi(p)
        else:
            def generator(idx):
                return a_series.coeff({
                    "alpha": idx,
                    "beta": i
                })
            self.a[(p, i)] = FPS(
                generator,
                alpha
            )
            self.a[(p, i)].name = "a_" + str(i) + " @ " + str(p)
        return self.a[(p, i)]
    
    def chi_power_op_cpnm1(self, p, n):
        assert n > 0
        ret = self.chi(p)**(p**n - 1) * self.get_x(p**n - 1) * p**(n-1)
        ret.name = "chi^" + str(2*(p**n-1)) + " P_" + str(p) + "(c_" + str(p**n - 1) + ")"
        return ret

    def chi_power_op_xpnm1(self, p, n):
        # mod p^m, ..., v_{n-1}
        m = p**n
        ret = self.chi_power_op_cpnm1(p, n) * Rational(ComplexCobordism.v(m), m)
        for i in range(1, n):
            d = p**i
            ret += self.chi_power_op_cpnm1(p, n - i) * self.chi_power_op_xpnm1(p, i)**(m//d) \
                * Rational(ComplexCobordism.mu(m, d) * ComplexCobordism.v(m), (m // d) * ComplexCobordism.v(d))
        ret = self.chromatically_truncate(ret, p, n)
        ret.name = "chi^" + str(2*(m-1)) + " P_" + str(p) + "(v_" + str(n) + ")"
        return ret
    
    def chromatically_truncate(self, f, p, n):
        # kill (p, ..., v_{n-1})
        def generator(idx):
            coeff = f.coeff(idx)
            if n == 0:
                return coeff

            for i in range(1, n):
                try:
                    coeff = coeff.subs(self.get_x(p**i - 1), 0)
                except AttributeError:
                    pass
            
            try:
                coeff = trunc(coeff, p)
            except polys.polyerrors.ComputationFailed:
                coeff = 0
            
            
            return coeff
        
        return FPS(
            generator,
            f.vars
        )
    
    def as_parseable_obj(self):
        obj = super().as_parseable_obj()
        obj["chis"] = {}
        for p in self.chis:
            obj["chis"][str(p)] = self.chis[p].as_parseable_obj()
        
        obj["a_series"] = {}
        for p in self.a_series:
            obj["a_series"][str(p)] = self.a_series[p].as_parseable_obj()
        
        return obj
    
    def from_parseable_obj(self, obj):
        super().from_parseable_obj(obj)
        if "chis" in obj:
            for p in obj["chis"]:
                if int(p) not in self.chis:
                    self.chis[int(p)] = self.chi(int(p))
                self.chis[int(p)].from_parseable_obj(obj["chis"][p])
        if "a_series" in obj:
            for p in obj["a_series"]:
                if int(p) not in self.a_series:
                    self.a_series[int(p)] = self.get_a_series(int(p))
                self.a_series[int(p)].from_parseable_obj(obj["a_series"][p])


MU = ComplexCobordism()
try:
    t0 = time.time()
    print("Loading MU")
    MU.load("./MU2.json")
    t1 = time.time()
    print("MU loaded in " + str(t1 - t0) + "s")
except FileNotFoundError:
    pass
