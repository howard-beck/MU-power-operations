from sympy import *
from FPS import *
from fgls_auto import *
from MU_fgl import *

p = 3

a_syms = {}
z = symbols("z")
def h_inv_gen(i):
    if i not in a_syms:
        a_syms[i] = symbols("a" + str(i))
    return a_syms[i]

n = 1
N = 3

#h_inv_sym = 
t0 = time.time()
chi = MU.chi(p)
chi_trunc = MU.chromatically_truncate(chi, p, N)
#chi.calculate_up_to(27)

P_v1 = MU.chi_power_op_xpnm1(p, n)
P_v1.calculate_up_to(17)

red_p_series = MU.get_n_series(p).shift(-1)
trunc_red_p_series = MU.chromatically_truncate(red_p_series, p, N)

k = (chi_trunc**(p**n - 1)).shift(-(p**n - 1))
#print("K:")
#k.calculate_up_to(10)
#MU.chromatically_truncate(k, p, N).calculate_up_to(12)

#new_P_v1 = MU.chromatically_truncate(P_v1, p, n) + trunc_red_p_series.shift((p-1)*(p**n - 1))# * k#ALPHA**((p - 2) * (p**n - 1))
new_P_v1 = MU.chromatically_truncate(P_v1, p, N) - k*trunc_red_p_series#ALPHA**((p - 2) * (p**n - 1))
trunc_red_p_series.name = "<" + str(p) + "> mod (p, ..., v_" + str(n-1) + ")"
#print("reduced p-series:")
#trunc_red_p_series.calculate_up_to(10)
#MU.chromatically_truncate(new_P_v1, p, n).calculate_up_to(10)

shifted_chiPvn = new_P_v1.shift(-2*(p-1) * (p**n - 1))
#chi_red = MU.chromatically_truncate(chi, p, n)
shifted_chi_power = (chi_trunc**(2*(p**n - 1))).shift(-(p-1) * 2*(p**n - 1))
shifted_chi_power_inv = shifted_chi_power.mult_inv(const_inv = 1)

the_power_op = shifted_chiPvn * shifted_chi_power_inv
the_power_op.name = "P(v" + str(n) + ")"
the_power_op = MU.chromatically_truncate(the_power_op, p, N)
print("=====")
#MU.chromatically_truncate(new_P_v1, p, n).calculate_up_to(10)
the_power_op.calculate_up_to(12)


#MU.chromatically_truncate(chi, p, 1).calculate_up_to(9)
t1 = time.time()
print("Performed calculations in " + str(t1 - t0) + "s")

#a_poly = FPS.prod(MU.get_beta_plus_n_series(i) for i in range(1, p))
#a_poly.calculate_up_to(5)

#MU.get_a(2, 0).calculate_up_to(5)
#h = 


t0 = time.time()
MU.save()
t1 = time.time()
print("Saved in " + str(t1 - t0) + "s")