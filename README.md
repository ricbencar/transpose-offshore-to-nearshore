# OFFSHORE-TO-NEARSHORE WAVE TRANSFORMATION

## Abstract

This repository implements a deterministic offshore-to-nearshore wave-transformation model for time-series forcing data stored in CSV format. The program reads offshore significant wave height, mean wave direction, and peak period, propagates each sea state to a prescribed local depth, and writes both a transformed time series and a statistical report. The implemented workflow combines linear wave dispersion, directional screening relative to a user-supplied coastline orientation, refraction by Snell-type kinematic transformation, shoaling through linear group-velocity relations, and a Miche-type breaking limitation.

The code is designed as a practical engineering post-processing tool rather than as a spectral wave model. It does not solve two-dimensional phase-resolving hydrodynamics, diffraction, current-wave interaction, or energy dissipation over complex bathymetry. Instead, it applies a compact and computationally efficient sequence of one-point transformations to each record independently. That makes the tool suitable for rapid forcing transposition, screening studies, preprocessing, sensitivity testing, and preparation of nearshore design series where a simplified but transparent physical model is acceptable.

This README is written to document the **actual implemented algorithm** in `transpose.cpp`, including numerical details, directional conventions, filtering rules, and reporting logic.

---

## 1. Scope and engineering purpose

The program transforms offshore wave conditions to a nearshore reference point at constant depth `depth_d` using a single-record, non-spectral formulation. For each valid time step, the code computes:

- deep-water wavelength;
- local wavelength from the linear dispersion relation;
- dimensionless depth `kh`;
- signed offshore obliquity relative to the coastline;
- refracted local obliquity;
- transformed local mean wave direction;
- shoaling coefficient;
- refraction coefficient;
- Miche breaking height;
- final local significant wave height capped by breaking.

In addition to the transformed series, the program writes a report with descriptive statistics and annual maxima.

The implementation is especially useful when one needs a reproducible engineering transformation chain that is:

- fast enough for large hindcast or reanalysis time series;
- explicit in its assumptions;
- simple to compile as a standalone executable;
- easy to audit because every computed quantity is written to the output.

---

## 2. Program outputs

Running the executable produces two files in the working directory:

### 2.1 `output.csv`

This file contains one transformed row per retained timestamp with the following columns:

```csv
datetime,swh_offshore,mwd_offshore,pp1d,L0,L,kh,alpha_offshore,alpha_local,
swh_local,mwd_local,Ks,Kr,Hb
```

### 2.2 `report.txt`

This file contains:

- the exact command line used to run the executable;
- descriptive statistics for all reported variables;
- a note explaining the filtering applied to selected local variables;
- annual maxima of `swh_offshore` and `swh_local`;
- overall maxima across all years.

---

## 3. Input data model

The executable expects a comma-separated file with at least four columns:

```csv
datetime,swh,mwd,pp1d
```

Additional columns may exist and are ignored.

### 3.1 Required input variables

- `datetime`: timestamp string; the code uses it as text and extracts the year from the first four characters;
- `swh`: offshore significant wave height;
- `mwd`: offshore mean wave direction;
- `pp1d`: offshore peak period.

### 3.2 Implied units and conventions

The implementation assumes:

- `swh` in metres;
- `pp1d` in seconds;
- `mwd` in degrees clockwise from North;
- `coast_dir` in degrees clockwise from North;
- `depth_d` in metres.

The program does not perform unit conversion. Inputs must already be expressed in a consistent system.

### 3.3 Direction convention

The code is written for the standard metocean convention in which mean wave direction is the **direction from which waves come**, expressed clockwise from geographic North. If the source dataset uses a going-to convention or a different angular reference system, the user must convert the data before running the program.

---

## 4. Meaning of `coast_dir`

`coast_dir` is the **coastline azimuth**, not the seaward normal. In other words, it represents the orientation of the shoreline itself, measured clockwise from North.

Examples:

- `coast_dir = 0` or `180` means a coastline aligned North-South;
- `coast_dir = 90` or `270` means a coastline aligned East-West.

This distinction is critical because the code computes wave obliquity relative to the coastline line, then infers the appropriate crest orientation and nearshore turning from that quantity.

---

## 5. Physical and mathematical formulation

### 5.1 Governing assumptions

The implementation follows the following simplified assumptions:

1. Each time step is treated independently.
2. Linear wave theory is used for wavelength, celerity, group velocity, shoaling, and refraction.
3. Refraction is represented through a local Snell-type relation based on the ratio `C/C0`.
4. Breaking is imposed by an upper bound derived from a Miche-type criterion.
5. Waves identified as arriving from the land side are suppressed and assigned zero local height.
6. No diffraction, bottom friction, wind input, currents, nonlinear triad interactions, or spectral spreading are modelled.

These assumptions make the tool intentionally compact. They also define its limits of validity.

---

### 5.2 Deep-water wavelength

For a wave period `T = pp1d`, the deep-water wavelength is computed as

$$
L_0 = \frac{g T^2}{2\pi},
$$

where:

- `g = 9.80665 m/s^2`;
- `T` is the offshore peak period.

In code, non-positive periods return `0`.

---

### 5.3 Local wavelength from the linear dispersion relation

At the target depth `d = depth_d`, the local wavelength `L` is obtained from the implicit linear dispersion equation rewritten as

$$
L = L_0 \tanh\left(\frac{2\pi d}{L}\right).
$$

This is equivalent to the classical dispersion relation

$$
\omega^2 = g k \tanh(k d),
$$

with

$$
L = \frac{2\pi}{k},
\qquad
\omega = \frac{2\pi}{T}.
$$

#### 5.3.1 Initial guess used in the code

The implementation does **not** start Newton iteration from a generic constant guess. It first computes

$$
k_0 d = \frac{2\pi d}{L_0},
$$

then forms a custom estimate for `kh`:

$$
(kh)_{\text{approx}} = \frac{k_0 d}{\tanh\left(\left(\frac{6}{5}\right)^{k_0 d}\sqrt{k_0 d}\right)}.
$$

From that quantity it builds the initial wavelength estimate as

$$
L_{\text{init}} = \frac{2\pi d}{(kh)_{\text{approx}}}.
$$

If this estimate is unusable, the code falls back to an Eckart-type approximation:

$$
L_{\text{Eckart}} = L_0 \sqrt{\tanh\left((k_0 d)^{1.25}\right)}.
$$

It then applies the safeguards:

- `L <= L0`;
- `L > 0`;
- if necessary, an additional fallback

$$
L = L_0 \tanh\left(\sqrt{k_0 d}\right),
$$

followed by `L = L0` as a final fallback.

#### 5.3.2 Newton iteration

After initialization, the code solves

$$
F(L) = L - L_0 \tanh\left(\frac{2\pi d}{L}\right) = 0
$$

with Newton-Raphson iteration:

$$
L_{n+1} = L_n - \frac{F(L_n)}{F'(L_n)}.
$$

The derivative used in the code can be written as

$$
F'(L) = 1 + L_0 \frac{2\pi d}{L^2} \left[1 - \tanh^2\left(\frac{2\pi d}{L}\right)\right].
$$

The iteration parameters hard-coded in the source are:

- maximum iterations: `20`;
- convergence tolerance: `1e-12`.

If an updated iterate becomes non-physical (`L <= 0`), iteration is stopped and the current estimate is retained.

---

### 5.4 Wavenumber and dimensionless depth

Once `L` is known, the local wavenumber is

$$
k = \frac{2\pi}{L},
$$

and the local depth parameter is

$$
kh = k d.
$$

If `L` is not positive, the code sets `k = 0` and therefore `kh = 0`.

---

### 5.5 Offshore directional screening relative to the coastline

The program first normalizes offshore wave direction to the interval `[0, 360)`.

It then computes a relative wave direction with respect to the coastline azimuth, also wrapped to `[0, 360)`.

In implementation terms, waves are treated as arriving from the **land side** when

$$
0 < \text{relativeDir} < 180.
$$

For those records, the code writes zero for all locally transformed quantities. This is a hard screening rule, not a gradual attenuation.

This behaviour is central to the model. It means the executable does not attempt to rotate land-side wave directions into a physically meaningful nearshore solution. Instead, it assumes such cases do not contribute to the local sea state of interest.

---

### 5.6 Signed offshore obliquity `alpha_offshore`

The source code computes a **signed crest-to-coast angle** called `alpha_offshore`. The crest orientation is inferred from the wave direction and coastline orientation.

If `relativeDir < 180`, the code uses

$$
\text{crest} = \text{mwd} - 90^{\circ},
$$

otherwise it uses

$$
\text{crest} = \text{mwd} + 90^{\circ}.
$$

The crest angle is wrapped into `[0, 360)`, then the signed difference with the coastline is reduced into `[-180, 180)`.

In practice, this creates a signed obliquity whose sign is later preserved during refraction.

---

### 5.7 Refraction and local angle `alpha_local`

The code applies a Snell-type transformation in terms of the angle relative to the coastline. Using the linear-wave celerity ratio

$$
\frac{C}{C_0} = \frac{L/T}{L_0/T} = \frac{L}{L_0} = \tanh(kh),
$$

the local angle satisfies

$$
\sin(\alpha_{\text{local}}) = \sin(\alpha_{\text{offshore}}) \tanh(kh).
$$

Therefore,

$$
\alpha_{\text{local}} = \arcsin\left[\sin(\alpha_{\text{offshore}})\tanh(kh)\right].
$$

The internal value passed to `asin` is clamped to `[-1, 1]` to avoid floating-point domain errors.

---

### 5.8 Local mean wave direction `mwd_local`

After refraction, the local mean wave direction is reconstructed by applying the angular change between offshore and local obliquity:

$$
\text{mwd}_{\text{local}} = \text{mwd}_{\text{offshore}} - \left(\alpha_{\text{offshore}} - \alpha_{\text{local}}\right),
$$

and the result is then wrapped to `[0, 360)`.

This is an implementation-level directional reconstruction. It preserves the code's sign convention and should be documented exactly because other formulations often reconstruct the local direction through the shore-normal or through a different angular basis.

---

### 5.9 Shoaling coefficient `Ks`

The shoaling coefficient is based on linear energy-flux conservation. The code uses

$$
K_s = \sqrt{\frac{C_{g0}}{C_g}},
$$

with the well-known finite-depth group velocity factor

$$
n = \frac{1}{2}\left(1 + \frac{2kh}{\sinh(2kh)}\right).
$$

Since

$$
\frac{C}{C_0} = \tanh(kh),
$$

the implemented expression becomes

$$
K_s = \sqrt{\frac{1}{\tanh(kh)\, 2n}}.
$$

Equivalently,

$$
K_s = \left[\tanh(kh)\left(1 + \frac{2kh}{\sinh(2kh)}\right)\right]^{-1/2}.
$$

#### 5.9.1 Numerical safeguards for `Ks`

The code returns `1.0` instead of evaluating the above formula when any of the following occur:

- `k <= 0`;
- `depth <= 0`;
- `sinh(2kh)` is too close to zero;
- `tanh(kh)` is too close to zero;
- the denominator becomes non-positive.

This is a practical numerical safeguard. It prevents singular or undefined values at extreme shallow-water or degenerate states.

---

### 5.10 Refraction coefficient `Kr`

The refraction coefficient is computed as

$$
K_r = \sqrt{\frac{\cos(\alpha_{\text{offshore}})}{\cos(\alpha_{\text{local}})}}.
$$

The code only evaluates this expression if

- `cos(alpha_offshore) >= 0`, and
- `cos(alpha_local) > tolerance`.

Otherwise, it falls back to

$$
K_r = 1.
$$

This matters numerically because grazing incidence can lead to very small denominators and unstable amplification.

---

### 5.11 Miche-type breaking wave height `Hb`

The breaking-limited height is computed as

$$
H_b = 0.142 \, L \, \tanh(kh).
$$

This is the expression explicitly implemented in the code and labelled as Miche (1944).

---

### 5.12 Local significant wave height `swh_local`

The unconstrained transformed height is

$$
H_{\text{trans}} = H_{s,\text{offshore}} K_s K_r.
$$

The final local significant wave height is then

$$
H_{s,\text{local}} = \min\left(H_{\text{trans}}, H_b\right).
$$

If `Hb <= 0`, the code uses `H_trans` directly. Any negative result is finally clipped to zero.

---

## 6. Algorithmic workflow

The executable follows the sequence below.

### 6.1 Input parsing

- Open input CSV.
- Read and discard the header line.
- Read all remaining non-empty lines.
- Ignore purely blank or whitespace-only rows.

### 6.2 Temporal ordering and duplicate removal

The program sorts records lexicographically by the first comma-separated field, interpreted as the datetime string. After sorting, duplicate timestamps are removed by retaining a single line per distinct datetime.

This has two practical consequences:

1. the output is sorted by datetime string, not by original file order;
2. if duplicate timestamps exist, only one record survives.

Because the sort is performed on the full list before deduplication and is not declared stable, duplicate rows with the same datetime should not be assumed to preserve original ordering.

### 6.3 Per-record transformation

Each retained record is parsed into:

- `datetime`;
- `swh`;
- `mwd`;
- `pp1d`.

If parsing fails, or if `swh <= 0`, or if `pp1d <= 0`, the code writes the original offshore values followed by zeros for all derived quantities.

For valid records, it then:

1. normalizes offshore direction;
2. checks whether the wave comes from the land side;
3. if land-side, sets local outputs to zero;
4. otherwise computes `L0`, `L`, `k`, `kh`, `alpha_offshore`, `alpha_local`, `mwd_local`, `Ks`, `Kr`, `Hb`, and `swh_local`.

### 6.4 Parallel execution

The record-processing loop is parallelized with OpenMP:

```cpp
#pragma omp parallel for schedule(static)
```

Thread-local containers are used for annual maxima, which are merged after the parallel loop. Output rows are first stored in memory and then written sequentially so that final file order matches the sorted order of retained timestamps.

---

## 7. Output variables and exact meanings

### 7.1 `swh_offshore`

The offshore significant wave height read from input.

### 7.2 `mwd_offshore`

The offshore mean wave direction wrapped to `[0, 360)`.

### 7.3 `pp1d`

The offshore peak period read from input.

### 7.4 `L0`

Deep-water wavelength from linear theory.

### 7.5 `L`

Finite-depth wavelength obtained from the Newton solution of the dispersion equation.

### 7.6 `kh`

Dimensionless depth parameter.

### 7.7 `alpha_offshore`

Signed offshore obliquity relative to the coastline.

### 7.8 `alpha_local`

Signed local obliquity after refraction.

### 7.9 `swh_local`

Final nearshore significant wave height after shoaling, refraction, and breaking limitation.

### 7.10 `mwd_local`

Local mean wave direction reconstructed from the obliquity change.

### 7.11 `Ks`

Shoaling coefficient.

### 7.12 `Kr`

Refraction coefficient.

### 7.13 `Hb`

Breaking-limited wave height.

---

## 8. Statistical report methodology

The generated `report.txt` is not a generic summary. It follows several precise rules implemented in the source.

### 8.1 Variables included in the report

The report covers the following thirteen variables:

1. `swh_offshore`
2. `mwd_offshore`
3. `pp1d`
4. `L0`
5. `L`
6. `kh`
7. `alpha_offshore`
8. `alpha_local`
9. `swh_local`
10. `mwd_local`
11. `Ks`
12. `Kr`
13. `Hb`

### 8.2 Linear statistics

For non-directional variables, the code computes:

- count;
- mean;
- sample standard deviation;
- minimum;
- percentiles at 1%, 10%, 25%, 50%, 75%, 90%, and 99%;
- maximum.

The percentile routine sorts the sample and performs linear interpolation between neighbouring ranks.

### 8.3 Hybrid circular statistics for directions

For `mwd_offshore` and `mwd_local`, the code uses a hybrid strategy:

- all angles are first wrapped to `[0, 360)`;
- circular mean is computed from the vector sum of unit phasors;
- circular standard deviation is computed from the mean resultant length;
- minimum, maximum, median, and all listed quantiles are computed linearly on the wrapped angles.

Thus, the mean and standard deviation are circular, while quantiles are not. This is an intentional design choice in the implementation and should be understood before interpreting the report.

### 8.4 Filtering applied to local variables

The following variables are **filtered** before descriptive statistics are computed:

- `alpha_local`
- `swh_local`
- `mwd_local`
- `Ks`
- `Kr`
- `Hb`

For these variables, the code excludes all records for which

$$
\text{swh}_{\text{local}} \le \text{tolerance}.
$$

This removes land-side waves and any time step whose final local wave height is effectively zero.

All other variables are summarized over the full retained record set.

### 8.5 Annual maxima

The code extracts the year as the first four characters of `datetime` and stores the annual maximum of:

- `swh_offshore`;
- `swh_local`.

The report also prints the overall maximum across all years.

---

## 9. Edge cases and implementation safeguards

The code contains several defensive choices that are important for correct interpretation.

### 9.1 Invalid or non-physical forcing

If any of the following apply:

- parsing failure;
- `swh <= 0`;
- `pp1d <= 0`;

then all derived quantities are written as zero.

### 9.2 Land-side waves

If

$$
0 < \text{relativeDir} < 180,
$$

the wave is treated as arriving from the land side and local outputs are set to zero.

### 9.3 Small-denominator protections

The code avoids unstable evaluations by checking for near-zero denominators in the shoaling and refraction calculations.

### 9.4 Non-physical Newton updates

If Newton iteration proposes `L <= 0`, the iteration stops immediately and the current estimate is retained.

### 9.5 Angle wrapping

All principal directional outputs are normalized either to `[0, 360)` or `[-180, 180)` depending on physical meaning:

- directions such as `mwd_offshore` and `mwd_local` are wrapped to `[0, 360)`;
- obliquities such as `alpha_offshore` are wrapped to `[-180, 180)`.

---

## 10. Practical interpretation of the model

This executable should be understood as a **point-transformation model**. It is well suited to cases where the user wants to map offshore forcing to a local depth under simplified assumptions and where the main priorities are transparency, reproducibility, and speed.

It is less suitable when the engineering problem requires any of the following:

- explicit bathymetric ray tracing over spatially varying depth;
- diffraction around breakwaters, headlands, or structures;
- current-induced refraction or Doppler shift;
- spectral partitioning or directional spreading;
- bottom-friction dissipation;
- surf-zone energy dissipation beyond a simple breaking cap;
- harbour resonance or basin-scale wave penetration.

In those cases, a spectral model such as SWAN or a more detailed nearshore propagation framework is more appropriate.

---

## 11. Usage

```sh
./transpose input_csv coast_dir depth_d
```

### 11.1 Arguments

- `input_csv`: input CSV file;
- `coast_dir`: coastline azimuth in degrees clockwise from North;
- `depth_d`: target local depth in metres.

### 11.2 Example

```sh
./transpose hindcast.csv 270 8.0
```

This example means:

- read `hindcast.csv`;
- use an East-West coastline orientation (`270` is equivalent to `90` as a line azimuth);
- transform waves to a local point at depth `8.0 m`.

---

## 12. Build and compilation

The source comment block documents the following compilation command:

```sh
g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ -o transpose transpose.cpp -lm
```

### 12.1 Meaning of the main flags

- `-O3`: high compiler optimization;
- `-fopenmp`: enables OpenMP parallelization;
- `-march=native`: tunes code generation to the compiling machine;
- `-std=c++17`: uses the C++17 language standard;
- `-Wall -Wextra -pedantic`: enables strict warning sets;
- `-Wconversion -Wsign-conversion`: requests warnings on implicit numeric conversions;
- `-static -static-libgcc -static-libstdc++`: requests static linkage where supported;
- `-lm`: links the math library.

---

## 13. Worked conceptual sequence

For one valid offshore record, the program conceptually applies the following chain:

$$
(T, H_s, \text{MWD})
\rightarrow
L_0
\rightarrow
L
\rightarrow
kh
\rightarrow
\alpha_{\text{offshore}}
\rightarrow
\alpha_{\text{local}}
\rightarrow
K_s, K_r
\rightarrow
H_b
\rightarrow
H_{s,\text{local}}, \text{MWD}_{\text{local}}.
$$

With equations grouped compactly, the implemented model is:

$$
L_0 = \frac{gT^2}{2\pi},
$$

$$
L = L_0 \tanh\left(\frac{2\pi d}{L}\right),
\qquad
k = \frac{2\pi}{L},
\qquad
kh = kd,
$$

$$
\sin(\alpha_{\text{local}}) = \sin(\alpha_{\text{offshore}})\tanh(kh),
$$

$$
K_s = \left[\tanh(kh)\left(1 + \frac{2kh}{\sinh(2kh)}\right)\right]^{-1/2},
\qquad
K_r = \sqrt{\frac{\cos(\alpha_{\text{offshore}})}{\cos(\alpha_{\text{local}})}},
$$

$$
H_b = 0.142 \, L \, \tanh(kh),
\qquad
H_{s,\text{local}} = \min\left(H_{s,\text{offshore}} K_s K_r, H_b\right).
$$

That compact system captures the entire implemented nearshore transformation for sea-side waves.

---

## 14. Limitations and recommended use

The model should be used with engineering judgment. In particular:

- it is best interpreted as a rapid transformation tool for a fixed target depth;
- it is appropriate for screening, long time-series preprocessing, and transparent design support;
- it should not be mistaken for a full coastal wave-propagation model;
- coastline orientation and directional convention must be checked carefully before use;
- the land-side suppression rule can strongly affect results and should be verified against local geometry;
- the use of peak period alone is a simplification when the incident sea state is broadband.

For important design decisions, results should be cross-checked against field knowledge, spectral modelling, or site-specific engineering assessment.

---

## 15. References underlying the implemented theory

The code structure reflects standard linear wave-transformation concepts commonly associated with:

- Airy wave theory for dispersion and celerity;
- Snell-type refraction of wave rays;
- linear group-velocity shoaling relations;
- Miche-type depth-limited breaking;
- Eckart-type finite-depth wavelength approximation for initialization.

This repository documents the implemented formulas directly from the source code, which should remain the primary reference for exact computational behaviour.

---

## 16. Repository summary

In operational terms, `transpose.cpp` is a compact engineering executable that:

- reads offshore wave forcing from CSV;
- removes duplicate timestamps after datetime sorting;
- transforms each record to a prescribed depth;
- rejects land-side waves;
- applies refraction, shoaling, and breaking limitation;
- writes a complete transformed series;
- produces descriptive statistics and annual maxima.

For users who need a transparent and auditable offshore-to-nearshore post-processing tool, this combination of physical clarity and implementation simplicity is the central value of the repository.