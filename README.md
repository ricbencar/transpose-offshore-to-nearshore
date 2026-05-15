# OFFSHORE-TO-NEARSHORE WAVE TRANSFORMATION

## Abstract

This repository implements a deterministic offshore-to-nearshore wave-transformation model for time-series forcing data stored in CSV format. The program reads offshore significant wave height, mean wave period, and mean wave direction, propagates each valid sea state to a prescribed local depth, and writes both a transformed time series and a statistical report.

The implemented workflow combines linear wave dispersion, directional screening relative to a user-supplied coastline orientation, refraction by a Snell-type kinematic transformation, shoaling through linear group-velocity relations, and a Miche-type breaking limitation.

The code is intended as a practical engineering post-processing tool rather than as a spectral wave model. It does not solve two-dimensional phase-resolving hydrodynamics, diffraction, wave-current interaction, spectral spreading, bottom friction, or energy dissipation over complex bathymetry. Instead, it applies a compact one-point transformation to each record independently.

This README documents the implemented behaviour of `transpose.cpp`, including input requirements, numerical formulation, directional conventions, filtering rules, output files, report logic, and compilation.

---

## 1. Scope and engineering purpose

The program transforms offshore wave conditions to a nearshore reference point at constant depth `depth_d`. The transformation is performed independently for each timestamp.

| Aspect | Implemented behaviour |
|---|---|
| Model type | One-point, deterministic, non-spectral transformation model |
| Input sea-state variables | `swh`, `mwp`, `mwd` |
| Required time variable | `datetime` |
| Geometry input | Coastline azimuth `coast_dir` and target depth `depth_d` |
| Main calculations | Linear dispersion, refraction, shoaling, breaking limitation |
| Directional convention | Wave direction is the direction **from which** waves come |
| Land-side waves | Suppressed by setting local transformed quantities to zero |
| Main outputs | `output.csv` and `report.txt` |

For each valid sea-side record, the code computes deep-water wavelength, finite-depth wavelength, dimensionless depth, offshore obliquity, local refracted obliquity, local mean wave direction, shoaling coefficient, refraction coefficient, Miche-type breaking height, and breaking-limited local significant wave height.

The tool is suitable for screening studies, sensitivity testing, preprocessing of offshore forcing, and preparation of simplified nearshore design series where transparency and reproducibility are more important than full spatial wave modelling.

---

## 2. Input CSV format

The executable expects a comma-separated file whose header starts with the following four columns:

```csv
datetime,swh,mwp,mwd
```

The program reads only these first four columns. Additional columns may be present after `mwd`, but they are ignored by the calculation and by the report. For example, the following extended header is valid:

```csv
datetime,swh,mwp,mwd,wind,dwi,u10,v10
```

| Column | Required position | Unit | Meaning | Used in calculations |
|---|---:|---|---|---|
| `datetime` | 1 | text | Timestamp string. The first four characters are used as the year for annual maxima. | Yes |
| `swh` | 2 | m | Offshore significant wave height. | Yes |
| `mwp` | 3 | s | Offshore mean wave period. This is the only period variable used by the wave transformation. | Yes |
| `mwd` | 4 | degrees | Offshore mean wave direction, clockwise from North. | Yes |
| Additional columns | 5 onward | any | Optional columns such as `wind`, `dwi`, `u10`, or `v10`. | No |

Header validation is case-insensitive and trims surrounding blanks or quotation marks. However, the first four field names must still resolve to `datetime`, `swh`, `mwp`, and `mwd` in this exact order.

The CSV parser is simple and comma-based. It is suitable for ordinary numeric CSV files but is not intended for fields containing embedded commas inside quoted text.

---

## 3. Units and directional convention

The program does not perform unit conversion. Inputs must already be expressed in a consistent system.

| Quantity | Expected convention |
|---|---|
| `swh` | metres |
| `mwp` | seconds |
| `mwd` | degrees clockwise from geographic North |
| `coast_dir` | degrees clockwise from geographic North |
| `depth_d` | metres |

Mean wave direction is interpreted using the standard metocean convention: it is the direction **from which** waves come. If the source dataset uses a going-to direction, or any other angular convention, the direction data must be converted before running the executable.

---

## 4. Output files

Running the executable creates two files in the working directory.

| File | Content |
|---|---|
| `output.csv` | Transformed time series, with one row per retained timestamp. |
| `report.txt` | Command line, descriptive statistics, local-variable filtering note, annual maxima, and overall maxima. |

The output CSV header is:

```csv
datetime,swh_offshore,mwd_offshore,mwp,L0,L,kh,alpha_offshore,alpha_local,swh_local,mwd_local,Ks,Kr,Hb
```

Numeric values in `output.csv` are written in fixed decimal notation with ten decimal places.

| Output variable | Meaning |
|---|---|
| `datetime` | Timestamp copied from the input row. |
| `swh_offshore` | Offshore significant wave height read from input column `swh`. |
| `mwd_offshore` | Offshore mean wave direction read from input column `mwd`, wrapped to `[0, 360)` for valid parsed records. |
| `mwp` | Offshore mean wave period read from input column `mwp`. |
| `L0` | Deep-water wavelength computed from `mwp`. |
| `L` | Finite-depth wavelength computed from `mwp` and `depth_d`. |
| `kh` | Dimensionless depth parameter. |
| `alpha_offshore` | Signed offshore crest-to-coast obliquity, in degrees. |
| `alpha_local` | Signed local obliquity after refraction, in degrees. |
| `swh_local` | Final local significant wave height after shoaling, refraction, and breaking limitation. |
| `mwd_local` | Local mean wave direction reconstructed after refraction. |
| `Ks` | Shoaling coefficient. |
| `Kr` | Refraction coefficient. |
| `Hb` | Miche-type breaking height. |

---

## 5. Command-line usage

The executable is called with exactly three arguments after the program name:

```sh
./transpose input_csv coast_dir depth_d
```

| Argument | Meaning |
|---|---|
| `input_csv` | Input CSV file whose first four columns are `datetime,swh,mwp,mwd`. |
| `coast_dir` | Coastline azimuth in degrees clockwise from North. |
| `depth_d` | Target local depth in metres. |

Example:

```sh
./transpose input.csv 270 8.0
```

This reads `input.csv`, uses a coastline azimuth of `270°`, and transforms the offshore records to a local depth of `8.0 m`.

The program stops with an error if `coast_dir` or `depth_d` cannot be parsed as numbers, or if `depth_d <= 0`.

---

## 6. Coastline orientation

`coast_dir` is the coastline azimuth, not the seaward normal. It represents the orientation of the shoreline itself, measured clockwise from North.

| `coast_dir` value | Coastline alignment |
|---:|---|
| `0` or `180` | North-South coastline |
| `90` or `270` | East-West coastline |

The code normalizes `coast_dir` to `[0, 360)` before using it. The distinction between coastline azimuth and shore-normal direction is important because the code computes wave obliquity relative to the coastline line and then reconstructs the local mean wave direction from the obliquity change induced by refraction.

---

## 7. Physical and mathematical formulation

The implementation is based on the following modelling assumptions.

| Assumption | Implementation consequence |
|---|---|
| Independent records | Each timestamp is transformed separately. |
| Linear wave theory | Used for wavelength, celerity, group velocity, shoaling, and refraction. |
| Snell-type refraction | Local angle is computed from the ratio between finite-depth and deep-water celerity. |
| Miche-type breaking cap | Local wave height is limited by a breaking height. |
| Land-side suppression | Waves identified as arriving from land are assigned zero local transformed quantities. |
| No spatial wave model | Diffraction, bathymetric ray tracing, currents, spectral spreading, bottom friction, and nonlinear interactions are not represented. |

### Deep-water wavelength

For a wave period `T = mwp`, the deep-water wavelength is:

$$
L_0 = \frac{gT^2}{2\pi},
$$

where `g = 9.80665 m/s²`. If `T <= 0`, the code returns `L0 = 0`.

### Finite-depth wavelength

At target depth `d = depth_d`, the finite-depth wavelength `L` is obtained from:

$$
L = L_0 \tanh\left(\frac{2\pi d}{L}\right).
$$

This is equivalent to the linear dispersion relation:

$$
\omega^2 = gk\tanh(kd),
\qquad
L = \frac{2\pi}{k},
\qquad
\omega = \frac{2\pi}{T}.
$$

The code first computes:

$$
k_0d = \frac{2\pi d}{L_0},
$$

then forms the initial estimate:

$$
(kh)_a = \frac{k_0d}{\tanh\left[\left(\frac{6}{5}\right)^{k_0d}\sqrt{k_0d}\right]}.
$$

From this value, the initial wavelength is:

$$
L_{\mathrm{init}} = \frac{2\pi d}{(kh)_a}.
$$

If this estimate is unusable, the code applies an Eckart-type fallback:

$$
L_{\mathrm{Eckart}} = L_0\sqrt{\tanh\left[(k_0d)^{1.25}\right]}.
$$

The estimate is capped so that `L <= L0`. If it remains unusable, the code falls back to:

$$
L = L_0\tanh\left(\sqrt{k_0d}\right),
$$

and finally to `L = L0` if needed.

Newton-Raphson iteration is then applied to:

$$
F(L) = L - L_0\tanh\left(\frac{2\pi d}{L}\right) = 0,
$$

with update:

$$
L_{n+1} = L_n - \frac{F(L_n)}{F'(L_n)},
$$

and derivative:

$$
F'(L) = 1 + L_0\frac{2\pi d}{L^2}\left[1 - \tanh^2\left(\frac{2\pi d}{L}\right)\right].
$$

The hard-coded iteration parameters are `20` maximum iterations and a convergence tolerance of `1e-12`. If a Newton update produces `L <= 0`, the iteration stops and the current estimate is retained.

### Wavenumber and dimensionless depth

Once `L` is known, the local wavenumber and dimensionless depth are:

$$
k = \frac{2\pi}{L},
\qquad
kh = kd.
$$

If `L` is not positive, the code sets `k = 0` and `kh = 0`.

### Directional screening relative to the coastline

The offshore mean wave direction is first normalized to `[0, 360)`. The code then computes the relative wave direction with respect to the coastline azimuth using the implemented expression:

```text
relativeDir = (mwd_offshore - coast_dir) mod 360
```

A wave is treated as arriving from the land side when:

```text
0 < relativeDir < 180
```

For those records, the code writes zero for all locally transformed quantities. This is a hard screening rule, not a gradual attenuation. The cases `relativeDir = 0` and `relativeDir >= 180` are treated as sea-side cases by the implementation.

### Signed offshore obliquity

The code computes a signed crest-to-coast angle called `alpha_offshore`. The wave crest orientation is inferred from the offshore mean wave direction as follows:

| Condition | Crest orientation used by the code |
|---|---|
| `relativeDir < 180` | `mwd - 90°` |
| `relativeDir >= 180` | `mwd + 90°` |

The crest angle is wrapped to `[0, 360)`. The signed difference between the crest and the coastline azimuth is then reduced to `[-180, 180)`. The sign of this obliquity is preserved during refraction.

### Refraction and local direction

The code applies a Snell-type relation using the linear-wave celerity ratio:

$$
\frac{C}{C_0} = \frac{L/T}{L_0/T} = \frac{L}{L_0} = \tanh(kh).
$$

Therefore:

$$
\sin(\alpha_{\mathrm{local}}) =
\sin(\alpha_{\mathrm{offshore}})\tanh(kh),
$$

and:

$$
\alpha_{\mathrm{local}} =
\arcsin\left[\sin(\alpha_{\mathrm{offshore}})\tanh(kh)\right].
$$

The internal value passed to the inverse sine is clamped to `[-1, 1]` to avoid floating-point domain errors.

The local mean wave direction is reconstructed using:

```text
mwd_local = wrap_to_360(mwd_offshore - (alpha_offshore - alpha_local))
```

This reconstruction preserves the sign convention used by the implemented obliquity calculation.

### Shoaling, refraction coefficient, and breaking cap

The shoaling coefficient is derived from linear energy-flux conservation:

$$
K_s = \left[\tanh(kh)\left(1 + \frac{2kh}{\sinh(2kh)}\right)\right]^{-1/2}.
$$

The refraction coefficient is:

$$
K_r = \sqrt{\frac{\cos(\alpha_{\mathrm{offshore}})}{\cos(\alpha_{\mathrm{local}})}}.
$$

The code only evaluates this expression when the offshore cosine is non-negative and the local cosine is greater than the numerical tolerance. Otherwise, `Kr = 1` is used as a safe fallback.

The Miche-type breaking height is:

$$
H_b = 0.142L\tanh(kh).
$$

The transformed height before breaking is:

$$
H_{\mathrm{trans}} = H_{\mathrm{offshore}}K_sK_r.
$$

The final local significant wave height is:

$$
H_{\mathrm{local}} = \min(H_{\mathrm{trans}}, H_b).
$$

If `Hb <= 0`, the code uses the transformed height directly. Any negative local height is clipped to zero.

---

## 8. Processing workflow

The program applies the following workflow.

| Stage | Implemented operation |
|---|---|
| Argument validation | Requires exactly `input_csv`, `coast_dir`, and `depth_d`; checks that numeric arguments are valid and that depth is positive. |
| Header validation | Requires the first four CSV header fields to be `datetime,swh,mwp,mwd`, after trimming and lower-casing. |
| Input parsing | Reads all non-empty rows after the header and parses the first four comma-separated fields. |
| Sorting | Sorts records lexicographically by the `datetime` field. |
| Duplicate removal | Removes duplicate timestamps after sorting, retaining one line per distinct timestamp. |
| Direction handling | Wraps valid offshore mean wave directions to `[0, 360)`. |
| Land-side screening | Sets local transformed quantities to zero when `0 < relativeDir < 180`. |
| Sea-side transformation | Computes `L0`, `L`, `kh`, `alpha_offshore`, `alpha_local`, `swh_local`, `mwd_local`, `Ks`, `Kr`, and `Hb`. |
| Parallel execution | Uses OpenMP for the per-record loop with `#pragma omp parallel for schedule(static)`. |
| Output writing | Stores output lines in memory and writes them sequentially so that final order follows the sorted and deduplicated input. |
| Annual maxima | Uses thread-local maps during parallel processing and merges them after the loop. |

Because sorting is performed on the full record list before duplicate removal and is not declared stable, duplicate rows with the same timestamp should not be assumed to preserve original file order.

---

## 9. Statistical report methodology

The report covers the following variables:

| Group | Variables |
|---|---|
| Offshore input variables | `swh_offshore`, `mwd_offshore`, `mwp` |
| Wave-theory variables | `L0`, `L`, `kh` |
| Directional transformation variables | `alpha_offshore`, `alpha_local`, `mwd_local` |
| Height-transformation variables | `swh_local`, `Ks`, `Kr`, `Hb` |

For non-directional variables, the code reports count, mean, sample standard deviation, minimum, 1st percentile, 10th percentile, 25th percentile, median, 75th percentile, 90th percentile, 99th percentile, and maximum. Percentiles are computed after sorting the sample and applying linear interpolation between neighbouring ranks.

For `mwd_offshore` and `mwd_local`, the code uses hybrid directional statistics: the mean and standard deviation are circular, based on unit-vector summation, while minimum, maximum, median, and percentiles are computed linearly on wrapped angles in `[0, 360)`.

The following local variables are filtered before descriptive statistics are computed:

```text
alpha_local, swh_local, mwd_local, Ks, Kr, Hb
```

For these variables, records are excluded from descriptive statistics when:

```text
swh_local <= 1e-12
```

This removes land-side waves and records whose final local wave height is effectively zero from the local-variable statistics. All other variables are summarized over the full retained record set.

Annual maxima are computed for `swh_offshore` and `swh_local`. The year is extracted from the first four characters of `datetime`. The report also prints an overall maximum row across all reported years.

---

## 10. Edge cases and safeguards

| Case | Behaviour |
|---|---|
| Wrong number of arguments | Program stops with a usage message. |
| Non-numeric `coast_dir` or `depth_d` | Program stops with an error. |
| `depth_d <= 0` | Program stops with an error. |
| Input file cannot be opened | Program stops with an error. |
| Output file cannot be created | Program stops with an error. |
| Missing header | Program stops with an error. |
| Header does not start with `datetime,swh,mwp,mwd` | Program stops with an error. |
| Valid header but no data rows | Program creates `output.csv` with only the header and creates `report.txt` stating that no valid data were present. |
| Row cannot be parsed | Derived variables are written as zero. |
| `swh <= 0` or `mwp <= 0` | Derived variables are written as zero. |
| Land-side wave | Local transformed quantities are written as zero. |
| Near-zero shoaling denominator | `Ks = 1` is used as a numerical safeguard. |
| Invalid refraction denominator | `Kr = 1` is used as a numerical safeguard. |
| Newton update gives `L <= 0` | Iteration stops and the current estimate is retained. |
| Direction output | Mean wave directions are wrapped to `[0, 360)`. |
| Obliquity output | Signed obliquities are wrapped to `[-180, 180)`. |

---

## 11. Build and compilation

The source is standard C++17 and uses OpenMP for parallel processing. A typical compilation command is:

```sh
g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ -o transpose transpose.cpp -lm
```

| Flag | Purpose |
|---|---|
| `-O3` | Enables high compiler optimisation. |
| `-fopenmp` | Enables OpenMP multithreading. |
| `-march=native` | Tunes the executable for the compiling machine. |
| `-std=c++17` | Compiles as C++17. |
| `-Wall -Wextra -pedantic` | Enables strict warning checks. |
| `-Wconversion -Wsign-conversion` | Reports implicit numeric-conversion warnings. |
| `-static -static-libgcc -static-libstdc++` | Requests static linkage where supported. |
| `-lm` | Links the mathematical library where required. |

On Windows with MSYS2/UCRT64, the same command can be used from the UCRT64 shell if `g++` and OpenMP support are installed.

---

## 12. Practical interpretation and limitations

This executable is a compact point-transformation model. It is appropriate when the objective is to transform a long offshore time series to a fixed local depth using a transparent, reproducible, and computationally light method.

| Suitable use | Not represented by this tool |
|---|---|
| Rapid offshore-to-nearshore transposition | Spatial bathymetric ray tracing |
| Long time-series preprocessing | Diffraction around structures or headlands |
| Screening and sensitivity tests | Wave-current interaction |
| Transparent design-support calculations | Directional spreading and spectral partitioning |
| Fixed-depth local transformation | Bottom-friction dissipation |
| Reproducible engineering post-processing | Surf-zone dissipation beyond the simple breaking cap |
| Auditable simplified calculations | Harbour resonance or basin-scale wave penetration |

For design-critical applications, results should be checked against engineering judgement, local bathymetric context, and, where appropriate, a spectral wave model or site-specific wave-transformation study.

---

## 13. Compact computational sequence

For one valid sea-side offshore record, the implemented calculation follows this sequence:

$$
(mwp, H_s, MWD)
\rightarrow
L_0
\rightarrow
L
\rightarrow
kh
\rightarrow
\alpha_{\mathrm{offshore}}
\rightarrow
\alpha_{\mathrm{local}}
\rightarrow
K_s, K_r
\rightarrow
H_b
\rightarrow
H_{\mathrm{local}}, MWD_{\mathrm{local}}.
$$

The main equations are:

$$
L_0 = \frac{gT^2}{2\pi},
\qquad
L = L_0\tanh\left(\frac{2\pi d}{L}\right),
\qquad
k = \frac{2\pi}{L},
\qquad
kh = kd,
$$

$$
\sin(\alpha_{\mathrm{local}}) = \sin(\alpha_{\mathrm{offshore}})\tanh(kh),
$$

$$
K_s = \left[\tanh(kh)\left(1 + \frac{2kh}{\sinh(2kh)}\right)\right]^{-1/2},
\qquad
K_r = \sqrt{\frac{\cos(\alpha_{\mathrm{offshore}})}{\cos(\alpha_{\mathrm{local}})}},
$$

$$
H_b = 0.142L\tanh(kh),
\qquad
H_{\mathrm{local}} = \min(H_{\mathrm{offshore}}K_sK_r, H_b).
$$

In this sequence, `T` is the input `mwp`.

---

## 14. Repository summary

In operational terms, `transpose.cpp` reads offshore wave forcing from CSV, requires the first four input columns to be `datetime,swh,mwp,mwd`, ignores any additional input columns, sorts records by timestamp, removes duplicate timestamps, transforms each valid sea-side record to a prescribed depth, suppresses land-side waves, applies refraction, shoaling, and breaking limitation, writes a complete transformed time series to `output.csv`, and writes descriptive statistics plus annual maxima to `report.txt`.

The program is therefore a compact engineering tool for transparent offshore-to-nearshore wave transposition based on mean wave period, significant wave height, mean wave direction, coastline azimuth, and target local depth.
