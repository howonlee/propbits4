import numpy as np
import numpy.fft as np_fft
import numpy.random as npr
import matplotlib.pyplot as plt
import functools
import time
import comhop
import impl

def a_kronecker_expansion_chart():
    num_propbits = 5
    rands = [npr.rand(2) for _ in range(num_propbits)]
    plt.clf()
    plt.title("Full Kronecker expansion of random propbitset")
    plt.plot(functools.reduce(np.kron, rands))
    plt.savefig("figs/a", bbox_inches="tight", dpi=300)

def inner_product(fst_pbset, snd_pbset):
    res = 1
    for fst_pb, snd_pb in zip(fst_pbset, snd_pbset):
        res *= np.dot(fst_pb, snd_pb)
    return res

def b_fast_fft_chart():
    ct_res, our_res = [], []
    propbit_sizes = range(10, 25)
    time.time()
    for num_propbits in propbit_sizes:
        pbs = [npr.randn(2) for _ in range(num_propbits)]
        before_our_res = time.time()
        fft_pbs = comhop.fft_propbits(num_propbits, 1)
        inner_product(fft_pbs, pbs)
        after_our_res = time.time()
        our_res.append(after_our_res - before_our_res)

    for num_propbits in propbit_sizes:
        pbs = [npr.randn(2) for _ in range(num_propbits)]
        kron_res = functools.reduce(np.kron, pbs)
        before_ct_res = time.time()
        curr_ct_res = np_fft.fft(kron_res)
        after_ct_res = time.time()
        ct_res.append(after_ct_res - before_ct_res)
    plt.clf()
    plt.title("One-Amplitude FFT Timings, \n Our method vs. Numpy method (standard Cooley-Tukey)")
    plt.semilogy(ct_res, label="Cooley-Tukey")
    plt.semilogy(our_res, label="Ours")
    plt.ylabel("Time (sec), log scale")
    plt.xlabel("Number propbits (log n)")
    plt.xticks(ticks=range(15), labels=range(10, 25))
    plt.legend(loc="upper left")
    plt.savefig("figs/b", bbox_inches="tight", dpi=300)


def c_complex_number_ansatz_for_modular_functions_chart():
    mod = 39
    mod_fn = (2 ** np.linspace(0, 99, 100)) % mod
    complex_fn = np.exp(2j * np.pi * (2 ** np.linspace(0, 99, 100) % mod) / mod)
    plt.clf()
    plt.title("Complex number ansatz for (2 ** x) % 39, x domain")
    plt.plot(mod_fn, label="Modulo")
    plt.plot(complex_fn, label="Complex")
    plt.legend(loc="upper left")
    plt.savefig("figs/c1", bbox_inches="tight", dpi=300)

    plt.clf()
    plt.title("Complex number ansatz for (2 ** x) % 39, frequency domain, adjusted scales")
    plt.plot(np.abs(np_fft.fft(mod_fn)[1:]) / 15, label="Modulo")
    plt.plot(np.abs(np_fft.fft(complex_fn)[1:]), label="Complex")
    plt.legend(loc="upper left")
    plt.savefig("figs/c2", bbox_inches="tight", dpi=300)

def d_trivial_spectrum_chart():
    mod = 39
    propbits = [[1, np.exp(2j * np.pi * (2 ** i) / mod)] for i in range(9, -1, -1)]
    propbits_fn = impl.prod_propbits(propbits)
    complex_fn = impl.full_fn_expand(propbits_fn, 10)

    fft_complex_fn = np_fft.fft(complex_fn)
    plt.clf()
    plt.title("Time space for trivial example, real part, x mod 39")
    plt.plot(complex_fn, label="Time value")
    plt.legend(loc="upper left")
    plt.savefig("figs/d1", bbox_inches="tight", dpi=300)

    plt.clf()
    plt.title("Frequency space for trivial example, normalized, x mod 39")
    plt.plot(np.abs(fft_complex_fn), label="Frequency value")
    plt.legend(loc="upper left")
    plt.savefig("figs/d2", bbox_inches="tight", dpi=300)

def e_vaguely_less_trivial_spectrum_chart():
    mod = 15
    propbits = [[1, np.exp(2j * np.pi * (2 ** i) / mod)] for i in range(9, -1, -1)]
    propbits[8] = [1, np.exp(2j * np.pi * (2 ** (10 - 8))/ (mod * 3))]
    propbits_fn = impl.prod_propbits(propbits)
    complex_fn = impl.full_fn_expand(propbits_fn, 10)

    fft_complex_fn = np_fft.fft(complex_fn)
    plt.clf()
    plt.title("Time space for marginally less trivial example,\nReal part, x mod 15")
    plt.plot(complex_fn, label="Time value")
    plt.legend(loc="upper left")
    plt.savefig("figs/e1", bbox_inches="tight", dpi=300)

    plt.clf()
    plt.title("Frequency space for marginally less trivial example,\nNormalized, x mod 15")
    plt.plot(np.abs(fft_complex_fn), label="Frequency value")
    plt.legend(loc="upper left")
    plt.savefig("figs/e2", bbox_inches="tight", dpi=300)

def f_entangled_spectrum_chart():
    mod1 = 39
    mod2 = 44
    propbits_1 = [[1, np.exp(2j * np.pi * (2 ** i) / mod1)] for i in range(9, -1, -1)]
    propbits_fn_1 = impl.prod_propbits(propbits_1)
    propbits_2 = [[1, np.exp(2j * np.pi * (2 ** i) / mod2)] for i in range(9, -1, -1)]
    propbits_fn_2 = impl.prod_propbits(propbits_2)
    complex_fn = impl.full_fn_expand(propbits_fn_1, 10) + impl.full_fn_expand(propbits_fn_2, 10)
    fft_complex_fn = np_fft.fft(complex_fn)
    plt.clf()
    plt.title("Time space for low-yaem entangled example,\nReal part")
    plt.plot(complex_fn, label="Time value")
    plt.legend(loc="upper left")
    plt.savefig("figs/f1", bbox_inches="tight", dpi=300)

    plt.clf()
    plt.title("Frequency space for low-yaem entangled example,\nNormalized")
    plt.plot(np.abs(fft_complex_fn), label="Frequency value")
    plt.legend(loc="upper left")
    plt.savefig("figs/f2", bbox_inches="tight", dpi=300)

    
def complex_modexp_repr(num_pbs, mod):
    # Do the below modulo mod if num_pbs is nontrivial.
    # Replace base if different base
    sequence = [pow(2, 2 ** n, mod) for n in range(num_pbs - 1, -1, -1)]

    fst = [np.exp(2j * np.pi / mod), np.exp(2j * np.pi * sequence[0] / mod)]
    rest = [[np.e, np.exp(sequence_member)] for sequence_member in sequence[1:]]
    res = [fst] + rest
    return np.array(res)

def g_complex_comhop_modexp_chart():
    mod = 51
    num_propbits = 10
    num_amplitudes = int(2 ** num_propbits)
    ordinary_modexp = np.array([pow(2, x, mod) for x in range(num_amplitudes)])
    complex_modexp = np.array([np.exp(2j * np.pi * pow(2, x, mod) / mod) for x in range(num_amplitudes)])
    comhop_pbs = complex_modexp_repr(num_propbits, mod)
    comhop_idx_fn = comhop.comhop_propbits(comhop_pbs)
    comhop_pb_res = impl.full_fn_expand(comhop_idx_fn, num_propbits)
    plt.clf()
    plt.title("Ordinary modular exponent versus complex ansatz, (2 ** x) mod 51, x domain")
    plt.plot(ordinary_modexp / 50, label="Modular exponent, scale adjusted", alpha=0.5)
    plt.plot(comhop_pb_res, label="Complex ansatz", alpha=0.5)
    plt.legend(loc="upper left")
    plt.savefig("figs/g1", bbox_inches="tight", dpi=300)
    plt.clf()
    plt.title("Ordinary modular exponent versus complex ansatz, (2 ** x) mod 51,\nOrdinary (ordinary ring) frequency domain")
    plt.plot(np.abs(np_fft.fft(ordinary_modexp))[1:] / 50, label="Modular exponent, scale adjusted", alpha=0.5)
    plt.plot(np.abs(np_fft.fft(comhop_pb_res))[1:], label="Complex ansatz", alpha=0.5)
    plt.legend(loc="upper left")
    plt.savefig("figs/g2", bbox_inches="tight", dpi=300)

def comhop_fft_apply(arr):
    fft_mat = np_fft.fft(np.eye(arr.shape[0]))
    return np.prod(np.exp(arr * fft_mat), axis=1)

def h_comhop_fft_coincidence_with_normal_fft_chart():
    # some complex crap, basically
    # comhop fft of that
    # fft of that 
    to_find = np.exp(2j * np.pi * (np.linspace(0, 499, 500) ** 2) / 21)
    normal_fft = np.abs(np_fft.fft(to_find)[1:])
    our_fft = np.abs(np.log(comhop_fft_apply(to_find)[1:]))
    plt.clf()
    plt.title("Coincidence of comhop-3 log FFT absval and ordinary FFT absval of \n exp(2j * pi * (x ** 2) / 21)")
    plt.plot(normal_fft, label="Ordinary FFT", alpha=0.5)
    plt.plot(our_fft, label="Comhop-3 log FFT", alpha=0.5)
    plt.legend(loc="upper left")
    plt.savefig("figs/h", bbox_inches="tight", dpi=300)

def i_comhop_mod_sq_chart():
    mod = 51
    mod_answer = (np.linspace(0, 499, 500) ** 2) % mod
    complex_answer = np.exp(2j * np.pi * (np.linspace(0, 499, 500) ** 2) / mod)
    plt.clf()
    plt.title("Coincidence of modular square and complex exponential square, x domain, scaled")
    plt.plot(mod_answer, label="Modulo", alpha=0.5)
    plt.plot((complex_answer * 25) + 25, label="Complex", alpha=0.5)
    plt.legend(loc="upper left")
    plt.savefig("figs/i1", bbox_inches="tight", dpi=300)
    
    plt.clf()
    plt.title("Coincidence of modular square and complex exponential square, frequency domain, scaled")
    plt.plot(np.abs(np_fft.fft(mod_answer)[1:]), label="Modulo", alpha=0.5)
    plt.plot(np.abs(np_fft.fft(complex_answer)[1:] * 25), label="Complex", alpha=0.5)
    plt.legend(loc="upper left")
    plt.savefig("figs/i2", bbox_inches="tight", dpi=300)

if __name__ == "__main__":
    a_kronecker_expansion_chart()
    b_fast_fft_chart()
    c_complex_number_ansatz_for_modular_functions_chart()
    d_trivial_spectrum_chart()
    e_vaguely_less_trivial_spectrum_chart()
    f_entangled_spectrum_chart()
    g_complex_comhop_modexp_chart()
    h_comhop_fft_coincidence_with_normal_fft_chart()
    i_comhop_mod_sq_chart()
