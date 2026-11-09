"""
Spin dynamics in a magnetic field with quadratic anisotropy.

H = -gamma * B . S + K_perp * Sx^2 + K_par * Sz^2      (hbar = 1)

Works for any spin s (1/2, 1, 3/2, ...). Basis ordering: m = s, s-1, ..., -s.
Time evolution: rho(t) = exp(-iHt) rho(0) exp(iHt),  <O>(t) = Tr[rho(t) O].

Usage:  python spin_dynamics.py            (runs s = 1/2, 1, 3/2)
        python spin_dynamics.py 1          (runs only s = 1)
"""
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# ------------------------------ Parameters ----------------------------------
gamma = 1.0
B = 6.0                      # field magnitude
K_perp = B * gamma * np.pi / 10
K_par = -0.1 * B * gamma
theta, phi = np.pi / 4, 0.0  # field direction (spin-1/2 uses theta = 0, see below)
t_max, n_points = 100.0, 20000   # time in ps


# ------------------------------ Spin operators ------------------------------
def spin_operators(s):
    """Sx, Sy, Sz for spin s (units of hbar), basis m = s, ..., -s."""
    m = np.arange(s, -s - 1, -1)
    dim = len(m)
    Sz = np.diag(m).astype(complex)
    Sp = np.zeros((dim, dim), dtype=complex)          # raising operator
    for k in range(1, dim):
        mk = m[k]
        Sp[k - 1, k] = np.sqrt(s * (s + 1) - mk * (mk + 1))
    Sm = Sp.conj().T
    Sx = (Sp + Sm) / 2
    Sy = (Sp - Sm) / (2j)
    return Sx, Sy, Sz


def hamiltonian(S, B, theta, phi):
    Sx, Sy, Sz = S
    H = -gamma * B * (np.sin(theta) * np.cos(phi) * Sx
                      + np.sin(theta) * np.sin(phi) * Sy
                      + np.cos(theta) * Sz)
    H = H + K_perp * Sx @ Sx + K_par * Sz @ Sz
    assert np.allclose(H, H.conj().T), "H is not Hermitian"
    return H


# ------------------------------ Initial states ------------------------------
def initial_states(s):
    """Dictionary name -> state vector (basis m = s, ..., -s)."""
    dim = int(round(2 * s + 1))
    states = {}
    for k in range(dim):
        v = np.zeros(dim, dtype=complex)
        v[k] = 1
        states[f"m = {s - k:g}"] = v
    if s == 0.5:
        states["+x direction"] = np.array([1, 1], dtype=complex) / np.sqrt(2)
    else:  # Schroedinger cat state (|s> + |-s>)/sqrt(2)
        cat = np.zeros(dim, dtype=complex)
        cat[0] = cat[-1] = 1 / np.sqrt(2)
        states["cat state"] = cat
    return states


# ------------------------------ Time evolution ------------------------------
def evolve(H, S, psi0, times):
    """Return <Sx>, <Sy>, <Sz>, |<S>| for all times."""
    E, V = np.linalg.eigh(H)
    c = V.conj().T @ psi0                         # coefficients in eigenbasis
    psi_t = V @ (c[:, None] * np.exp(-1j * E[:, None] * times[None, :]))
    exp = [np.einsum('it,ij,jt->t', psi_t.conj(), O, psi_t).real for O in S]
    mag = np.sqrt(exp[0]**2 + exp[1]**2 + exp[2]**2)
    return exp[0], exp[1], exp[2], mag


# ------------------------------ Plotting ------------------------------------
def run(s):
    S = spin_operators(s)
    # Spin-1/2: field along z only (as in the report); anisotropy is trivial for s = 1/2
    th = 0.0 if s == 0.5 else theta
    H = hamiltonian(S, B, th, phi)
    times = np.linspace(0, t_max if s != 0.5 else 10, n_points)
    labels = [r'$\langle S_x\rangle$', r'$\langle S_y\rangle$',
              r'$\langle S_z\rangle$', r'$|\langle \vec S\rangle|$']
    for name, psi0 in initial_states(s).items():
        data = evolve(H, S, psi0, times)
        fig, axes = plt.subplots(4, 1, figsize=(12, 9), sharex=True)
        for ax, y, lab in zip(axes, data, labels):
            ax.plot(times, y, lw=1.2)
            ax.set_ylabel(lab, fontsize=13)
            ax.grid(True, alpha=0.4)
        axes[-1].set_xlabel('Time (ps)', fontsize=12)
        fig.suptitle(f'Spin-{s:g}, initial state: {name}  '
                     rf'($B={B:g}$, $\theta={th/np.pi:.2f}\pi$, $\phi={phi/np.pi:.2f}\pi$)',
                     fontsize=13)
        plt.tight_layout()
        fname = f"figures/spin{s:g}_{name.replace(' ', '').replace('=', '').replace('/', '_')}.png"
        plt.savefig(fname, dpi=120)
        plt.close(fig)
        print("saved", fname)


if __name__ == "__main__":
    import os
    os.makedirs("figures", exist_ok=True)
    spins = [float(a) for a in sys.argv[1:]] or [0.5, 1, 1.5]
    for s in spins:
        run(s)
