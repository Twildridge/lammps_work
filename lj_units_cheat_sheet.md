# LJ Units Cheat Sheet for Swollen Polymer Membranes

---

## 1. Fundamental LJ Scales

All quantities in `units lj` derive from three independent scales:

| Scale | Symbol | Definition | PEG/water mapping |
|-------|--------|------------|-------------------|
| Length | σ | Bead diameter | 0.76 nm (PEG Kuhn length) |
| Energy | ε | LJ well depth | 4.14 × 10⁻²¹ J (k_B T at 300 K) |
| Mass | m | Bead mass | 1.54 × 10⁻²⁵ kg (93 g/mol, PEG Kuhn segment) |

**Convention:** T* = k_B T / ε = 1.0 means the simulation temperature equals the
physical temperature used to define ε. This is the standard Kremer–Grest choice.

---

## 2. Derived Units

| Quantity | LJ expression | Physical value (PEG) |
|----------|---------------|----------------------|
| Time | τ = σ √(m/ε) | 4.6 ps |
| Pressure | P* = ε/σ³ | 9.4 MPa |
| Force | F* = ε/σ | 5.4 pN |
| Diffusivity | D* = σ²/τ | 1.26 × 10⁻⁹ m²/s |
| Viscosity | η* = √(mε)/σ² | 1.68 × 10⁻⁴ Pa·s |

---

## 3. Simulation Parameters (Current Values)

### Pair interactions

| Pair | ε | σ | r_c | Style | Notes |
|------|---|---|-----|-------|-------|
| Default (all) | 1.0 | 1.0 | 2.5σ | Full LJ | Attractive tail |
| Polymer–solvent | εSP | 1.0 | 2.5σ | Full LJ | Tunable; currently 1.0 (θ) |
| Solvent–solvent | εSS | 1.0 | 2.5σ | Full LJ | Reference |
| Polymer–support/piston | 1.0 | 1.0 | 1.22σ | WCA | Purely repulsive |
| Solvent–support | 0 | — | 0 | Off | Solvent passes through support |
| Solvent–piston | 1.0 | 1.0 | 1.22σ | WCA | Active in equilibration/permeation |

### FENE bonds (Kremer–Grest standard)

```
bond_coeff * 30.0 1.5 1.0 1.0
             K    R₀   ε   σ
```

| Parameter | Value | Physical meaning |
|-----------|-------|------------------|
| K | 30 ε/σ² | Spring constant |
| R₀ | 1.5 σ | Maximum bond extension |
| Equilibrium bond length | 0.96 σ | Set by FENE + WCA balance |
| Max chain extensibility | 1.5σ per bond | Prevents unphysical stretching |

**Reference:** Kremer & Grest, J. Chem. Phys. 92, 5057 (1990).

### Thermostat and barostat

| Parameter | LJ value | Physical |
|-----------|----------|----------|
| T* | 1.0 | 300 K |
| P* (NPT target) | 0.1 | 0.9 MPa ≈ 9 atm |
| τ_T (thermostat damping) | 1.0 τ | 4.6 ps |
| τ_P (barostat damping) | 1.0 τ | 4.6 ps |
| dt (production) | 0.005 τ | 23 fs |

### Atom types

| Type | Species | Mass | Notes |
|------|---------|------|-------|
| 1 | Crosslinker | 1.0 | Diamond lattice nodes |
| 2 | Chain bead | 1.0 | Interior beads along chains |
| 3 | Solvent | 1.0 | Represents ~20 water molecules |
| 4 | Support | 1.0 | Frozen rigid body (bottom) |
| 5 | Piston | 1.0 | Rigid body, z-mobile |
| 6 | Walls | 1.0 | Frozen rigid body (lateral) |

---

## 4. Pressure–Density Relationship (LJ EOS at T* = 1.0)

NPT sets the pressure; the barostat adjusts box volume to find ρ*.

| P* | ρ* | Physical pressure | Physical ρ | Notes |
|-----|------|-------------------|------------|-------|
| 0.01 | 0.84 | 0.1 MPa (1 atm) | ~1.0 g/cm³ | Ambient; near vapor–liquid coexistence |
| **0.1** | **0.85** | **0.9 MPa (9 atm)** | **~1.0 g/cm³** | **Current setting** |
| 0.5 | 0.87 | 4.7 MPa | — | Moderate compression |
| 1.0 | 0.88 | 9.4 MPa (93 atm) | — | Previous setting (too high) |
| 5.0 | 0.98 | 47 MPa | — | Heavily compressed |

**Warning:** P* < ~0.01 at T* = 1.0 risks crossing into the vapor–liquid coexistence
region, causing the box to explode. P* = 0.1 is safely liquid-phase.

**Reference:** Johnson, Zollweg, Gubbins, Mol. Phys. 78, 591 (1993).

---

## 5. Piston Force → Applied Pressure

The piston is a dense hexagonal monolayer (spacing 0.2σ) spanning the box cross-section.

| Quantity | Expression | Value |
|----------|------------|-------|
| Piston face area | A ≈ L_x × L_y | ~87 × 87 ≈ 7600 σ² (post-NPT) |
| N_piston atoms | A / (spacing² × √3/2) | ~63,000 |
| Force per atom (start) | 0.01 ε/σ | 0.054 pN |
| Force per atom (end) | 0.07 ε/σ | 0.38 pN |
| Total force (end) | N × 0.07 | ~4400 ε/σ |
| Applied pressure (end) | F_total / A | ~0.58 ε/σ³ ≈ **5.5 MPa** |

### Reasonable force ranges

| Experiment type | Typical pressure | P* range | Force/atom range |
|-----------------|------------------|----------|------------------|
| Hydrogel mechanical testing | 0.1–10 MPa | 0.01–1.0 | 0.001–0.07 |
| Membrane permeation (RO) | 1–7 MPa | 0.1–0.7 | 0.01–0.07 |
| Gentle osmotic compression | 1–100 kPa | 0.0001–0.01 | 0.0001–0.001 |
| Confined compression (AFM) | 0.01–1 MPa | 0.001–0.1 | 0.0001–0.01 |

---

## 6. Solvent Quality and the χ Parameter

| εSP | Solvent quality | χ (mean-field, z ≈ 10) | Physical analog |
|-----|-----------------|------------------------|-----------------|
| < 0.95 | Poor | > +0.5 | Collapsed gel / LCST above transition |
| 0.95–0.99 | Marginal | +0.1 to +0.5 | PEG in water near ~40 °C |
| **1.0** | **θ (athermal)** | **0** | **Current setting; clean baseline** |
| 1.01–1.05 | Good | −0.1 to −0.5 | PEG in water at 25 °C |
| > 1.05 | Very good | < −0.5 | Strong swelling |

**Mean-field formula** (rough guide only; unreliable at φ_p > 0.1):

$$\chi \approx \frac{z}{2}\left(\varepsilon_{PP} + \varepsilon_{SS} - 2\varepsilon_{SP}\right), \quad z \approx 8\text{–}12$$

**Better approach:** Measure equilibrium φ_p at several εSP values, fit to Flory–Rehner
to extract effective χ.

### Why θ conditions (εSP = εSS = 1.0) for baseline

At θ, every solvent atom has the same per-atom virial regardless of whether its
neighbors are polymer or solvent. The LAMMPS partial stress then maps cleanly
to the continuum partial stress from mixture theory:

- σ_{s,zz}(inside gel) / σ_{s,zz}(outside) ≈ φ_s(inside) / φ_s(outside)
- Pore pressure p_p is spatially uniform at equilibrium (∇p_p = 0)
- The stress step at the gel boundary provides a clean reference for measuring
  deviations under compression/flow

When εSP ≠ εSS, the per-atom virial depends on local composition, creating an
offset between the LAMMPS bin-averaged quantity and the continuum pore pressure.

---

## 7. Mapping to Physical Hydrogel Systems

| Property | PEG hydrogel | Polyacrylamide | This simulation |
|----------|-------------|----------------|-----------------|
| Kuhn length (σ) | 0.76 nm | 0.56 nm | 1.0 (LJ) |
| Kuhn segment mass | 93 g/mol | 40 g/mol | 1.0 (LJ) |
| C∞ | 4.1 | 8.5 | — |
| Monomers per Kuhn segment | 2.1 | 1.1 | — |
| χ (water, 25 °C) | 0.40–0.45 | 0.48 | 0 (θ baseline) |
| Swollen φ_p | 0.01–0.3 | 0.02–0.2 | ~0.5 (measured) |
| Typical M (modulus) | 10–1000 kPa | 1–100 kPa | TBD |

**Note:** The high φ_p ≈ 0.5 in this simulation reflects the dense diamond lattice
(short chains, high crosslink density). This is more representative of a tightly
crosslinked membrane than a loosely swollen hydrogel. The modulus will be
correspondingly higher.

**Kuhn length references:**
Rubinstein & Colby, *Polymer Physics* (2003), Table 2.1;
Devanand & Selser, Macromolecules 24, 5943 (1991).

---

## 8. Quick Conversion Table

| To convert | Multiply LJ value by | To get |
|------------|---------------------|--------|
| Length (σ → nm) | 0.76 | nm |
| Energy (ε → J) | 4.14 × 10⁻²¹ | J |
| Energy (ε → kJ/mol) | 2.49 | kJ/mol |
| Time (τ → ps) | 4.6 | ps |
| Pressure (P* → MPa) | 9.4 | MPa |
| Pressure (P* → atm) | 93 | atm |
| Force (F* → pN) | 5.4 | pN |
| Temperature (T* → K) | 300 | K |
| Density (ρ* → g/cm³) | ~1.18 | g/cm³ (for PEG/water) |
| Timesteps (10⁶ at dt=0.005 → ns) | 23 | ns |

---

## 9. Key Conventions and Pitfalls

**FENE R₀:** Standard Kremer–Grest is R₀ = 1.5σ. Using R₀ = 2.5 allows
unphysical chain extensibility and was a previous error in this project.

**NPT pressure:** P* = 1 is ~93 atm, not ambient. Use P* = 0.1 for a safe
liquid-phase state near standard conditions. P* < 0.01 risks vapor formation.

**Mass equality:** All bead masses = 1.0. Fine for equilibrium properties
(modulus, stress profiles, φ). Incorrect for dynamics (D_coop, permeability).
Solvent beads should be ~4× heavier than polymer for PEG/water dynamics.

**Initial bond length:** The generator places beads at ~0.52σ (compressed).
A soft push-off phase (harmonic bonds + pair_style soft) is required before
switching to FENE with R₀ = 1.5 to avoid divergence.

**Stress computation:** LAMMPS `stress/atom` returns stress × volume in units
of pressure × volume (ε). Divide by bin volume to get pressure. The isotropic
partial stress is −(σ_xx + σ_yy + σ_zz)/(3 V_bin).

**Pore pressure definition** (from mixture theory):

$$p_p = \lambda + \frac{\partial A}{\partial \phi_s}$$

where λ is the Lagrange multiplier (hydraulic pressure) and ∂A/∂φ_s contains
osmotic/mixing contributions. At equilibrium: ∇p_p = 0 (uniform chemical potential).
In LAMMPS: p_p = σ_{s,zz} / φ_s (measured from binned partial stress and volume fraction).