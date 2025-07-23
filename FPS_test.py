from sympy import *
from FPS import *
from fgls_auto import *
from MU_fgl import *

p = 2

a_syms = {}
z = symbols("z")
def h_inv_gen(i):
    if i not in a_syms:
        a_syms[i] = symbols("a" + str(i))
    return a_syms[i]
    
#h_inv_sym = 
t0 = time.time()
chi = MU.chi(p)

n = 2
P_v1 = MU.chi_power_op_xpnm1(p, n)
P_v1.calculate_up_to(10)

red_p_series = MU.get_n_series(p).shift(-1)
trunc_red_p_series = MU.chromatically_truncate(red_p_series, p, n)

new_P_v1 = P_v1 + trunc_red_p_series * ALPHA**((p - 2) * (p**n - 1))
#MU.chromatically_truncate(new_P_v1, p, n).calculate_up_to(10)

PP_v1 = new_P_v1.shift(-2*(p**n - 1))
print("=====")
MU.chromatically_truncate(PP_v1, p, n).calculate_up_to(10)
#MU.chromatically_truncate(chi, p, 1).calculate_up_to(9)
t1 = time.time()
print(t1 - t0)

#a_poly = FPS.prod(MU.get_beta_plus_n_series(i) for i in range(1, p))
#a_poly.calculate_up_to(5)

#MU.get_a(2, 0).calculate_up_to(5)
#h = 


MU.save()