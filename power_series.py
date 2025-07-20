from sympy import *
from utils import *

alpha = symbols("alpha")
beta = symbols("beta")
x = symbols("x")

default_acc = 10



def get_accuracy(ps, var = alpha):
    acc = 1
    while ps.coeff(var, acc+1) != 0:
        acc += 1
    return acc

def monic_power_series_inverse(ps, var = alpha, acc = None):
    assert ps.coeff(var, 0) == 0 or ps.coeff(var, 0) == 0.0
    assert ps.coeff(var, 1) == 1 or ps.coeff(var, 1) == 1.0

    if acc is None:
        acc = get_accuracy(ps)
    # calculate the inverse of a power series with constant term 0 and linear term 1
    # using the Lagrange inversion theorem

    f = [ps.coeff(var, n) for n in range(0, acc+1)]

    g = [0] * (acc+1)
    g[1] = 1

    # get coefficients f[2], ...
    args = [f[i] for i in range(2, len(f))]
    # calculate derivatives
    derivs = [f[i] * factorial(i) for i in range(len(f))]
    f_hat = [Rational(1, i) * derivs[i] for i in range(2, len(f))]

    for n in range(2, acc + 1):
        total = 0
        for k in range(1, n):
            # calculate rising factorial n^{\underline{k}}
            rising = 1
            for i in range(0, k):
                rising *= (n + i)
            # Lagrange inversion theorem
            term = (-1)**k * rising * bell(n-1, k, f_hat)
            total += term
        
        # n-th coefficient of inverse
        g[n] = simplify(Rational(1, factorial(n)) * total)
    return sum(g[n] * var**n for n in range(acc+1))

def mult_inv(ps, var = alpha, const_inv = None):
    # get multiplicative inverse of a power series
    if const_inv is None:
        const_inv = Rational(1, a[0])
    a = [ps.coeff(var, i) for i in range(N+1)]
    b = [const_inv]

    for n in range(1, N+1):
        b.append(-const_inv(1, a[0]), sum(a[j] * b[n - j] for j in range(1, n+1)))
    
    return sum(b[i]*var**i for i in range(N+1))

def truncate(fs, var = alpha, M = default_acc):
    # truncate a power series up to N
    return sum([fs.coeff(var, i)*var**i for i in range(M+1)])

def double_truncate(fs, var1, var2, M1 = default_acc, total = False, M2 = None):
    if not total:
        if M2 is None:
            M2 = M1
    
    def get_bound(i):
        if total:
            return M1 - i
        else:
            return M2
        
    var1_truncs = [truncate(fs.coeff(var1, i), var2, get_bound(i)) for i in range(M1+1)]
    return sum([expand(var1_truncs[i] * var1**i) for i in range(M1+1)])

def mult_trunc(fs, var = alpha, M = default_acc):
    gs = [f + 0*var for f in fs]
    rets = [1 + 0*var, gs[0]]
    for i in range(1, len(fs)):
        prod = 0
        for deg in range(M+1):
            for j in range(deg+1):
                prod += expand(rets[-1].coeff(var, j) * gs[i].coeff(var, deg - j) * var**deg)
        rets.append(prod)
    return rets

def mult_trunc2(fs, var1, var2, M1 = default_acc, total = False, M2 = None):
    rets = [1, fs[0]]
    for i in range(1, len(fs)):
        print(str(i) + "/" + str(len(fs)))
        prod = expand(rets[-1] * fs[i])
        rets.append(double_truncate(prod, var1, var2, M1 = M1, total = total, M2 = M2))
        print(str(i) + "/" + str(len(fs) - 1))
    return rets

def subs_trunc(f, arr, var, M):
    g = expand(f + 0*var)

    taskbar = WIPBar(
        len(arr)
    )
    for (var_sub, val_sub) in arr:
        val_pows = mult_trunc([val_sub]*M, var, M)

        g_new = 0
        for i in range(M+1):
            c = g.coeff(var_sub, i)
            g_new += mult_trunc([c, val_pows[i]])[-1]
        g = g_new

        taskbar.update()
    return g