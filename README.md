# Molecular Dynamics of Triatomic Argon — NVT + SHAKE

_(for the subject: Mètodes Avançats de Simulació Molecular)_

_(by: itxasoma (Itxaso Muñoz Aldalur) and ManelDC55 (Manel Díaz Calvo))_

Classical MD simulation of 200 triatomic Argon molecules in the NVT ensemble.
Each molecule is a rigid equilateral triangle (side r₀ = 3 Å) enforced via the SHAKE algorithm.
Intermolecular interactions use a truncated Lennard-Jones potential (cutoff 2.5σ).
Temperature is controlled with a Berendsen thermostat.

## Structure


```
.
├── src/
│   ├── parameters.f90     # Global parameters (GLOBAL module)
│   ├── force.f90          # LJ force + PBC
│   ├── motion.f90         # Velocity Verlet integrator + SHAKE
│   ├── io_teacher.f90     # I/O routines
│   ├── rdf_module.f90     # Radial distribution function
│   ├── main.f90           # Program: free liquid (no constraints)
│   ├── mainshake.f90      # Program: constrained system (SHAKE)
│   ├── plots.ipynb        # Analysis and plots (requires python virtual environment)
│   ├── plots.gnu          # Testing plots in gnuplot
│   └── Makefile
│   ├── inputs/
│   │   ├── conf.data             # Initial configuration
│   │   ├── confnova.data         # Initial configuration
│   │   ├── exercicishake.dim     # Dimensions for exercicishake.f program
│   │   └── exercicishake.dades   # Simulation parameters
├── out/                   # Output files (created at runtime)
├── plots/                 # Plots (created by notebook)
├── mplstyle/              # Plotting style template
├── profe/
│   ├── exercicishake.f     # Main reference program, given by the professor, adapted as module program in src/
│   └──  gr.f90
└──           
```


## Compilation

```bash
cd src/       # enter the source file
make          # builds both executables
make clean    # removes object files and executables
```

## Running
```bash
make run        # free liquid (no SHAKE) → out/energies_T.dat, trajectory.xyz, gr_ArAr.dat
make runshake   # with SHAKE            → out/energies_T_shake.dat, trajectory_shake.xyz, gr_ArAr_shake.dat
```
Requires `gfortran`.
For the plots, `python3` or more. See requirements at first cell of the notebook.

## Parameters (`inputs/exercicishake.dades`)

| Parameter | Value | Units |
|-----------|-------|-------|
| N_mol | 200 | — |
| σ | 3.405 | Å |
| ε | 119.8 | K |
| m | 39.948 | g/mol |
| r₀ (SHAKE) | 3.0 | Å |
| T_ref | 225 | K |
| Δt | 0.02 | ps |
| τ_T | 0.1 | ps |
| N_steps | 10 000 | — |

## Output

| File | Description |
|------|-------------|
| `energies_T.dat` | t, Upot, K, Etot, T, λ (no SHAKE) |
| `energies_T_shake.dat` | same, with SHAKE |
| `trajectory.xyz` | atomic positions (no SHAKE) |
| `trajectory_shake.xyz` | atomic positions (SHAKE) |
| `gr_ArAr.dat` | g(r) intermolecular (no SHAKE) |
| `gr_ArAr_shake.dat` | g(r) intra+inter (SHAKE) |

## References

- Ryckaert, Ciccotti & Berendsen, *J. Comput. Phys.* **23**, 327 (1977) — SHAKE
- Berendsen et al., *J. Chem. Phys.* **81**, 3684 (1984) — thermostat
- Rahman, *Phys. Rev.* **136**, A405 (1964) — Ar LJ parameters
```



