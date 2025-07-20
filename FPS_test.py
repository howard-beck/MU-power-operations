from sympy import *
from FPS import *


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

use_hazewinkel = True

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
    if i not in xs:
        xs[i] = symbols("x" + str(i))
    return xs[i]

def get_v(i):
    if i not in vs:
        vs[i] = symbols("v" + str(i))
    return vs[i]

def get_c(i):
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

alpha = symbols("alpha")

def log_fgl_generator(n):
    if n == 0:
        return 0
    if n == 1:
        return 1
    return get_c(n-1)/n

log_fgl = FPS(
    log_fgl_generator,
    vars = alpha,
    name = "log"
)

exp_fgl = log_fgl.comp_inv()

for i in range(20):
    print(exp_fgl.get_coeff(i))