from sympy import *
import math
from power_series import *
from fgls import *



import random

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

xs = {}
vs = {}
cs = {}

accuracy = 20

use_hazewinkel = True
if use_hazewinkel:
    for i in range(1, accuracy+1):
        xs[i] = symbols("x" + str(i))
else:
    for i in range(1, accuracy+1):
        cs[i] = symbols("c" + str(i))

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
    vd = v(d)
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
            ret *= c(p, d)
    return ret

def get_x(i):
    if i == 0:
        return 1
    if i > accuracy:
        return 0
    if i not in xs:
        xs[i] = symbols("x" + str(i))
    return xs[i]

def get_v(i):
    if i > accuracy:
        return 0
    if i not in vs:
        vs[i] = symbols("v" + str(i))
    return vs[i]

def get_c(i):
    if i > accuracy:
        return 0
    if i not in cs:
        if i == 0:
            cs[0] = 1
        else:
            if use_hazewinkel:
                k = i + 1

                ret = Rational(k, v(k)) * get_x(k-1)

                for d in range(2, k):
                    if k % d == 0:
                        coeff = Rational(mu(k, d) * d, v(d))
                        ret += coeff * get_x(d-1)**(int(k/d+0.00001)) * get_c(int(k/d+0.00001)-1)
                
                cs[i] = ret
            else:
                cs[i] = symbols("c" + str(i))
    return cs[i]


#MU_fgl = FGL(
#    log_fgl = sum([get_c(i-1) * Rational(1, i) * alpha**i for i in range(1, accuracy+1)]), name = "MU"
#)
#MU_fgl.save()
try:
    MU_fgl = FGL.load("./MU.json")
except FileNotFoundError:
    MU_fgl = FGL(
        log_fgl = sum([get_c(i-1) * Rational(1, i) * alpha**i for i in range(1, accuracy+1)]), name = "MU"
    )
if MU_fgl.accuracy < accuracy:
    print("Not enough accuracy available, recalculating formal group law")
    MU_fgl = FGL(
        log_fgl = sum([get_c(i-1) * Rational(1, i) * alpha**i for i in range(1, accuracy+1)]), name = "MU"
    )
MU_fgl.save()

p = 2

p_series = MU_fgl.get_n_series(p)

chi = mult_trunc([MU_fgl.get_n_series(i) for i in range(1, p)])[-1]



x = symbols("x")


def coeffs_of_prod(fs, var, n):
    if len(fs) == 1:
        return [fs[0].coeff(var, i) for i in range(n+1)]
    else:
        coeffs_of_next_prod = coeffs_of_prod(fs[1:], var, n)
        ret = [
            sum(mult_trunc([fs[0].coeff(var, i), coeffs_of_next_prod[deg - i]], alpha, accuracy)[-1] for i in range(deg+1))
            for deg in range(n+1)
        ]
        return ret

x_plus_n_series_terms = [MU_fgl.get_x_plus_n_series(i) for i in range(p)]
ai = coeffs_of_prod(x_plus_n_series_terms, x, int(accuracy / 2))[1:]
ai_sym = [symbols("a" + str(i)) for i in range(len(ai))]

chi2 = symbols("chi2")
h_coeffs = [chi2**(-1)]

for n in range(1, len(ai)):
    h_coeffs.append(sum(mult_trunc([-ai_sym[j], chi2**(-1), h_coeffs[n - j]], alpha, accuracy)[-1] for j in range(1, n+1)))
    print(str(n) + "/" + str(len(ai)))

h = sum(h_coeffs[k] * x**k for k in range(len(h_coeffs)))
for i in range(len(h_coeffs)):
    print("h[x^" + str(i) + "]: " + str(h_coeffs[i]))

# now work mod p
# kill up to v_{HEIGHT}
HEIGHT = 0
def chromatic_reduction(poly):
    for n in range(1, HEIGHT + 1):
        print("Killing " + str(get_x(p**n - 1)))
        poly = poly.subs(get_x(p**n - 1), 0)
    if HEIGHT > -1:
        return trunc(poly, p)
    return poly

#h = chromatic_reduction(h)

print("Produced h")
def get_chi_power_op_cn(m):
    if int(accuracy/2 - 1 + 0.0001) < m:
        print("h is only known up to order " + str(int(accuracy/2 - 1)) + ", " + str(m) + " requested")
        print("Accuracy must be at least " + str(2*(m+1)))
        assert False
    h_to_mp1 = mult_trunc2([h]*(m+1), x, alpha, m, accuracy)[-1]

    chi_power_P_cm = chi2**(2*m+1) * sum(
        get_c(m - k) * h_to_mp1.coeff(x, k)
        for k in range(m+1)
    )

    ret = expand(chi_power_P_cm)

    return ret

def get_chi_power_op_xn(m):
    k = m + 1

    # v_{k-1}
    ret = get_chi_power_op_cn(m)

    for d in range(2, k):
        if k % d == 0:
            coeff = Rational(mu(k, d) * d, v(d))
            ret += coeff * mult_trunc([get_chi_power_op_xn(d-1)**(int(k/d+0.00001)), get_chi_power_op_cn(int(k/d+0.00001)-1)])[-1]
    
    ret *= Rational(v(m), m)

    return ret

def get_power_op_xn(m):
    assert p == 2
    chi_power_op = subs_trunc(
        get_chi_power_op_xn(m),
        [(ai_sym[i], ai[i]) for i in range(len(ai))] + [(chi2, chi)],
        var = alpha,
        M = accuracy
    )

    print(chi_power_op)
    print("")
    print(chromatic_reduction(chi_power_op))

get_power_op_xn(3)