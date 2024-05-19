
from dataclasses import dataclass
from functools import reduce
from math import sqrt, factorial
from typing import Any


@dataclass
class WmmModel:
    title: str
    date: str
    epoch: float
    nm: list[tuple[int, int]]
    g: list[float]
    h: list[float]
    g_dot: list[float]
    h_dot: list[float]


def diag_index(n: int, m: int) -> int:
    return n * (n + 1) // 2 + m - 1


def S_n_m(n: int, m: int) -> float:
    def double_fac(i: int) -> int:
        return reduce(lambda a, b: a * b, range(1, n + 1, 2))
    kron_m_0 = (1 if m == 0 else 0)
    return sqrt(
        ((2 - kron_m_0) * factorial(n - m)) / factorial(n + m)
    ) * double_fac(2 * n - 1) / factorial(n - m)


def print_number(x: float, max_width: int) -> str:
    return f"REAL({x:> {max_width}.16e}"

def gen_coeff_table(nm: list[tuple[int, int]], g: list[float], h: list[float]) -> str:
    max_width = max(max(len(str(x)) for x in col) for col in (g, h))
    def print_num(x: float): return print_number(x, max_width)
    code_lines: list[str] = ["// Auto-generated table by `tools/gen_coeffs.py`"]
    for (n_i, m_i), g_i, h_i in zip(nm, g, h):
        line = f"{{ .g = {print_num(g_i)}, .h = {print_num(h_i)} }},  // (n = {n_i:>3}, m = {m_i:>3})"
        code_lines.append(line)
    return "\n".join(code_lines)


def read_wmm_cof(fname: str) -> WmmModel:
    with open(fname) as f:
        # Parse header
        header = [x.strip() for x in f.readline().strip().split()]
        assert len(header) == 3
        epoch, name, date = header
        epoch = float(epoch)

        ns: list[int] = []
        ms: list[int] = []
        gs: list[float] = []
        hs: list[float] = []
        g_dots: list[float] = []
        h_dots: list[float] = []
        for line in f:
            line = line.strip()
            if all(x == "9" for x in line):
                break
            data = line.split()
            assert len(data) == 6

            # Order of the coefficients
            n, m = (int(x) for x in data[:2])
            ns.append(n)
            ms.append(m)

            # Compute correctly normed coefficients
            coeffs_i = (float(x) for x in data[2:])
            normed_coeffs = (S_n_m(n, m) * x for x in coeffs_i)
            normed_coeffs = (0.0 if x == 0.0 else x for x in normed_coeffs)
            g, h, g_dot, h_dot = normed_coeffs
            gs.append(g)
            hs.append(h)
            g_dots.append(g_dot)
            h_dots.append(h_dot)

    # Converts 2D coefficients (n, m) to 1D using a specific diagonalization index scheme
    N = len(ns)
    coeffs: Any = [None] * N
    for n, m, g, h, g_dot, h_dot in zip(ns, ms, gs, hs, g_dots, h_dots):
        ind = diag_index(n, m)
        coeffs[ind] = ((n, m), g, h, g_dot, h_dot)
    assert all(x is not None for x in coeffs)

    nm, g, h, g_dot, h_dot = zip(*coeffs)
    return WmmModel(name, date, epoch, nm, g, h, g_dot, h_dot)  # type: ignore
