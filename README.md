# Higher spin dynamics in the presence of anisotropy

Python simulation of the time evolution of spin-1/2, spin-1 and spin-3/2 systems in a magnetic field with quadratic anisotropy. Spin-1 and spin-3/2 are "qudits", which are of interest for quantum information. Project done at Uppsala University (specialized course topic in quantum technology, spring 2026), supervised by Erik Sjöqvist.

## Model

Hamiltonian (ħ = 1):

H = −γB·S + K⊥ Sx² + K∥ Sz²

The field direction is set by a polar angle θ and an azimuthal angle φ. The state evolves according to the Liouville–von Neumann equation, ρ(t) = e^(−iHt) ρ(0) e^(iHt). We compute ⟨Sx⟩, ⟨Sy⟩, ⟨Sz⟩ and |⟨S⟩| as functions of time by diagonalizing H numerically.

Parameters: K⊥/γB = π/10, K∥/γB = −0.1, θ = π/4, φ = 0.

## Main results

- **Spin-1/2** (test case): Larmor precession, with ⟨Sz⟩ = 0 and constant |⟨S⟩|, as expected analytically.
- **Spin-1**: clear periodic patterns (full quantum revivals) for the initial state m = 0 and for the Schrödinger cat state (|−1⟩ + |+1⟩)/√2. The states m = ±1 show no clear pattern.
- **Spin-3/2**: no clear revival for any tested initial state (m = ±1/2, ±3/2, cat state), which indicates interference between incommensurate phase factors.

Plots are in the `figures/` folder, and the full write-up is in `report_project_anisotropy(1).pdf`.



## Files

| File | Description |
|------|-------------|
| `spin_dynamics.py` | Simulation for any spin (1/2, 1, 3/2), all initial states |
| `figures/` | Generated plots |
| `report_project_anisotropy(1).pdf` | Project report |

## How to run

```
pip install numpy matplotlib
python spin_dynamics.py        # runs spin-1/2, 1 and 3/2
python spin_dynamics.py 1      # runs only spin-1
```

All initial states (each m state, plus the cat state for s ≥ 1 and the +x state for s = 1/2) are run automatically and saved in `figures/`. Field strength, anisotropy and field angles are set at the top of the script.

## Reference

Anticoherent states: D. Baguette and J. Martin, *Anticoherence measures for pure spin states*. Qudits: R. Bertlmann and P. Krammer, *Bloch vectors for qudits*.
