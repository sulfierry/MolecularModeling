# Molecular dynamics input parameters

This document outlines common parameters used in molecular dynamics simulations, covering stages of minimization, relaxation, and production.

## Common parameters in all stages

- `imin`: Indicates minimization execution. `imin=1` for energy minimization.
- `ntx`: Reads coordinates from the "inpcrd" file. Supported options: `1` (standard format without initial velocity) and `2`.
- `irest`: Restarts simulation. `irest=0` (default) does not restart.
- `ntpr`: Prints progress every `ntpr` steps.
- `ntwx`: Writes coordinates to the `mdcrd` file every `ntwx` steps. `ntwx=0` does not record.
- `cut`: Sets non-bonded cutoff in Angstroms. Common value: `8.0` for PME.

## Minimization

```
imin=1
ntx=1
irest=0
maxcyc=2000 # Maximum cycles of minimization
ncyc=1000 # Stages of steepest descent before conjugate gradient
ntpr=100
ntwx=0
cut=10.0
```

> `maxcyc`: Maximum cycles of minimization. `maxcyc=2000` for `xmin` minimizer with TNCG method.

### LMOD procedure

- Uses eigenvectors of low-frequency vibrational modes.
- Starts with energetically minimized molecular model.
- Utilizes ARPACK for vibrational modes.
- `Hv` calculated via finite differences: `Hv = [∇(xmin+h)−∇(xmin)]/h`.
- Reuses eigenvectors for efficiency in LMOD.

### XMIN
float xmin( float func(), int natm, float x[], float g[], float ene, float grms_out, struct xmod_opt xo);


- `xmin()`: Returns minimized energy and updates coordinates to the minimum energy conformation.

## General settings for relaxation, pre-Production, and production

```
imin=0
dt=0.002
ntf=2
ntc=2
temp0=300.0
ntpr=100
ntwx=2000
ntwr=2000
cut=10.0
ntb=2
ntp=1
ntt=2
nmropt=1
ig=-1
iwrap=1
ntr=1
restraintmask="X"
```


### Explained parameters

- `dt`: Time step. Maximum recommended: `0.002` with SHAKE.
- `ntf` and `ntc`: Force evaluation and bond length constraints with SHAKE.
- `tempi`: Initial temperature.
- `temp0`: Reference temperature. Standard: `300K`.
- `ntwx` and `ntwr`: Frequency of recording coordinates and restart file.
- `ntb`: Control of periodic boundaries. Standard for constant pressure: `ntb=2`.
- `ntp`: Constant pressure dynamics. `ntp=1` for isotropic position scaling.
- `ntt`: Temperature scaling. `ntt=2` for canonical ensemble. Note that setting `ntt=0` corresponds to the microcanonical ensemble `(NVE)` (which should approximate the canonical for large numbers of degrees of freedom). Some aspects of the "weak coupling ensemble" `(ntt=1)` have been examined and approximately interpolated between the microcanonical and canonical ensembles. Options `ntt=2` and `3` correspond to the canonical ensemble (constant T).
- `ntt=2` as we define the Andersen-type temperature coupling scheme, in which "collisions" randomize velocities to a distribution corresponding to `temp0` every `vrand` steps. Note that between these "massive collisions," the dynamics are Newtonian. Therefore, time correlation functions (etc.) may be calculated in these sections and the results calculated into an initial canonical distribution. Also, note that a very high collision rate (a very small value of `vrand`) will slow the rate at which molecules explore configuration space, while a very low rate means that the canonical distribution of energies will be sampled slowly. A discussion of this rate is provided by `Andersen`. Note that this option is not equivalent to the original thermostat described by `Andersen`.
- `nmropt`: Control of NMR parameters and variations. `nmropt=1` for constraints and changes.
- `ig`: Seed for pseudo-random number generator. `-1` for seed based on date/time.
- `iwrap`: Packing coordinates in the primary box. `iwrap=1` for packing. This means that for each molecule, its nearest periodic image to the middle of the "primary box" (with x coordinates between 0 and a, y coordinates between 0 and b, and z coordinates between 0 and c) will be the one written to the output file. This usually makes the resulting structures look visually better but has no effect on the energy or forces. Running such packing, however, can disrupt diffusion and other calculations.
- For very long runs, it may be necessary to set `iwrap = 1` to prevent coordinate output from overflowing the trajectory and restart file formats, especially if trajectories are written in the `ASCII` format instead of `NetCDF` (see also the `ioutfm` option).
- `ntr`: Cartesian space atom restraint using harmonic potential. `ntr=1` with `restraintmask` set.

---

This format provides a clear and concise structure, making it easier to understand the parameters and their functions in molecular dynamics simulations.


---

Esse formato proporciona uma estrutura clara e concisa, facilitando a compreensão dos parâmetros e suas funções nas simulações de dinâmica molecular.
