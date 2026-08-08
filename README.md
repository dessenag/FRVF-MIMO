# A Tutorial for MIMO Modal Identification via FRVF

**Authors:** Gabriele Dessena¹, Marco Civera², Beatrice E. Bauret Martínez¹

¹ Department of Aerospace Engineering, Universidad Carlos III de Madrid, Spain
² Department of Structural, Geotechnical and Building Engineering, Politecnico di Torino, Italy

This repo contains the code for the preprint:

B. E. Bauret Martínez, G. Dessena, M. Civera, and O. E. Bonilla-Manrique, "Enhanced input stacking for non-square MIMO modal identification of aeronautical structures via Fast and Relaxed Vector Fitting," [arXiv:2605.16037](https://arxiv.org/abs/2605.16037), 2026. Preprint.

When using this code for your work or research please cite the following:

1. B. E. Bauret Martínez, G. Dessena, M. Civera, and O. E. Bonilla-Manrique, "Enhanced input stacking for non-square MIMO modal identification of aeronautical structures via Fast and Relaxed Vector Fitting," [arXiv:2605.16037](https://arxiv.org/abs/2605.16037), 2026. Preprint.
2. B. E. Bauret Martínez, G. Dessena, M. Civera, and O. E. Bonilla-Manrique, "Multi-Input Multi-Output Fast and Relaxed Vector Fitting for Aircraft Ground Vibration Testing," Engineering Proceedings, vol. 133, no. 1, p. 162, 2026, doi: [10.3390/engproc2026133162](https://doi.org/10.3390/engproc2026133162).

**A note on attribution:** Fast and Relaxed Vector Fitting (FRVF) is not a new algorithm introduced by this work. The underlying Vector Fitting method, and the "fast" and "relaxed" refinements that give FRVF its name, are due to B. Gustavsen and co-authors — see the **Licence** section below for the required citations — and are implemented in `vectfit3.m`, third-party code from SINTEF Energy Research. What this repository and the reference article above contribute is the enhanced input-stacking strategy that extends Vector Fitting to non-square Multi-Input Multi-Output (MIMO) systems.

The tutorial illustrates, via the **noiseless numerical beam case of Section III** of [1], the capability of Fast and Relaxed Vector Fitting (FRVF) — extended with an enhanced input-stacking strategy for non-square Multi-Input Multi-Output (MIMO) systems — to extract modal parameters from Ground-Vibration-Test-like data, benchmarked against the classical Least-Squares Complex Exponential (LSCE) method and screened via stabilisation diagrams that automatically retain the highest stable model order. It does not reproduce the noise-robustness sweep or the experimental BAE Systems Hawk T1A aircraft case of the preprint; for those, see [1] and its associated data/software repositories.

This tutorial repository is structured after, and indebted to, [iLF-MIMO](https://github.com/dessenag/iLF-MIMO), the equivalent tutorial for the improved Loewner Framework (iLF).

## Repository Structure

| Path | Description |
|---|---|
| [`FRVF.m`](/FRVF.m) | Low-level Fast and Relaxed Vector Fitting engine (enhanced input stacking; wraps `vectfit3.m`) |
| [`FRVF_id.m`](/FRVF_id.m) | Modal-parameter extraction from MIMO FRVF |
| [`tutorial.mlx`](/tutorial.mlx) | MATLAB Live Script walking through the noiseless Section III beam case end-to-end, in MATLAB live-document format |
| [`tutorial.pdf`](/tutorial.pdf) | The above exported in a PDF document |
| [`tutorial_plain.m`](/tutorial_plain.m) | The same tutorial as a classic, traditionally-commented `.m` script — for version control, screen readers, terminals, and MATLAB releases without Live Editor |
| [`Utilities/lsce_fr.m`](/Utilities/lsce_fr.m) | Least-Squares Complex Exponential (LSCE) baseline identification (via MATLAB `modalfit`) |
| [`Utilities/stabilisation_diagram.m`](/Utilities/stabilisation_diagram.m) | Stabilisation-diagram construction and order-by-order stability screening; interactive pole-picking is optional and off by default — both tutorials take the highest stable order automatically instead |
| [`Utilities/compute_mac.m`](/Utilities/compute_mac.m) | Modal Assurance Criterion (MAC) between two mode-shape sets |
| [`Utilities/beam_element.m`](/Utilities/beam_element.m), [`assemble_beam.m`](/Utilities/assemble_beam.m) | 3-D Euler–Bernoulli beam finite-element assembly for the tutorial's example structure |
| [`Utilities/classify_modes.m`](/Utilities/classify_modes.m), [`mode_match_str.m`](/Utilities/mode_match_str.m) | Tutorial display/formatting helpers (strong-/weak-axis labelling, Table I row matching) |
| [`Utilities/get_vectfit3.m`](/Utilities/get_vectfit3.m) | Checks for the third-party `vectfit3.m` dependency and obtains it on first run (automatic download, with a manual fallback) — **see licence note below** |

## Requirements

MATLAB with the **Signal Processing Toolbox** (required by `lsce_fr.m` for `modalfit`, and by both tutorials for `pwelch`). Developed against recent MATLAB releases (R2021b or later recommended).

## Quick start

Open `tutorial.mlx` in MATLAB for the formatted live-document experience, or `tutorial_plain.m` for a plain script that runs identically. Running either from the repository root (or with the repository added to the path, which both do automatically via `addpath(genpath(...))`) reproduces Table I of [1] for the noiseless MIMO beam at the highest stable model order, the FRF-fit-quality plot at minimum order, and the stabilisation diagrams for both FRVF and LSCE.

On first run, either tutorial (and `FRVF.m`, if called directly) will check for the third-party `vectfit3.m` dependency; if it isn't found, you'll be walked through downloading it — see **Licence** below.

Both tutorials contain identical logic and call the same `Utilities/` functions (MATLAB's Live Editor does not support local function definitions in plain-text Live Code files, so the beam finite-element assembly and display helpers live there rather than at the bottom of either tutorial). `tutorial_plain.m` is complete and independently verified by running it end to end in MATLAB; `tutorial.mlx` is a genuine binary Live Script generated via MATLAB's own Editor (`matlab.desktop.editor.openDocument`/`saveAs`), independently checked by unzipping it and inspecting `matlab/document.xml` directly for correctly structured headings, bold text, equations, and the file/role table (rather than raw markup text).

## Licence

The original code in this repository (root `.m` files and `Utilities/`) is released under the **GNU General Public License v3.0** — see [LICENSE](/LICENSE).

`vectfit3.m` (required by `FRVF.m`) is **third-party code by Bjørn Gustavsen (SINTEF Energy Research)**, restricted to **non-commercial use only**, and is therefore **not bundled** in this repository — it is incompatible with the GPL-3.0 above. [`Utilities/get_vectfit3.m`](/Utilities/get_vectfit3.m) is called automatically by `FRVF.m` and both tutorials the first time `vectfit3.m` is needed. It always shows the licence restriction and the official SINTEF source page first, then asks you to explicitly accept before doing anything:

- If you accept, it downloads `VFIT3.zip` directly from SINTEF and extracts `vectfit3.m` into `Utilities/` itself, discarding the bundled example scripts/papers it doesn't need.
- If you decline, or the automatic download fails for any reason (network issue, moved/changed file, etc.), it falls back to a manual flow: full instructions are printed in the Command Window first, then the source page opens in your browser, then you're asked to locate the file you downloaded yourself.

`vectfit3.m`'s licence requires citing the following whenever it is used, reproduced here for convenience (see also the References section of [`Utilities/get_vectfit3.m`](/Utilities/get_vectfit3.m)):

1. B. Gustavsen and A. Semlyen, "Rational approximation of frequency domain responses by Vector Fitting," IEEE Trans. Power Delivery, vol. 14, no. 3, pp. 1052–1061, 1999.
2. B. Gustavsen, "Improving the pole relocating properties of vector fitting," IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587–1592, 2006.
3. D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter, "Macromodeling of Multiport Systems Using a Fast Implementation of the Vector Fitting Method," IEEE Microwave and Wireless Components Letters, vol. 18, no. 6, pp. 383–385, 2008.

## Known issues

`stabilisation_diagram.m` uses a Windows-style path separator when `saveFigures = true` (figure export to `ID_RESULTS\...`); on macOS/Linux, replace it with `fullfile('ID_RESULTS', ...)` before enabling exports. Both tutorials call the function with `saveFigures = false` throughout, so this does not affect the tutorials themselves.

## Repository Citation

B. E. Bauret Martínez, G. Dessena, M. Civera, and O. E. Bonilla-Manrique, "Enhanced input stacking for non-square MIMO modal identification of aeronautical structures via Fast and Relaxed Vector Fitting," [arXiv:2605.16037](https://arxiv.org/abs/2605.16037), 2026. Preprint.
