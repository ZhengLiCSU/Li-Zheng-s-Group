# Tunnel leakage reliability under freeze–thaw cycling

Companion data and code for peer review of *Reliability analysis of tunnel leakage under freeze–thaw cycling and spatial variability: Coupled thermo-hydro-mechanical modelling and hydraulically informed surrogates*.

[Code guide](CODE_GUIDE.md) describes the three stages: the physical solver, the hydraulic inputs and the fitted Kriging predictor.

## Reproduce the reported checks

From this directory, use Python with NumPy and SciPy. Choose a fresh output directory outside the checkout.

```sh
python -B verify_review.py --output-dir ../tunnel_leakage_review_results
python -B code/verify_threshold_mesh.py
```

The first command checks file integrity, reconstructs the 28 hydraulic and statistical inputs, compares the fitted Kriging predictions with the stored predictions for 80 held-out cases, and checks the Monte Carlo sample-size comparisons. It also regenerates the histogram data for Fig. 12(d). The second command re-integrates the accepted flow histories for the threshold-adjacent mesh comparison.

## Find the evidence

| Review question | Files |
| --- | --- |
| THM implementation and numerical controls | `solver/source/`, `solver/example/`, `solver/numerical_controls.json` |
| Random fields and Fig. 8 | `data/example/`, `code/plot_field_and_appendix.py` |
| Kriging implementation and Fig. 12 | `FEATURE_DEFINITIONS.txt`, `data/model/`, `data/validation/`, `data/figure12/`, `code/plot_figure12.py` |
| Mesh and Monte Carlo sensitivity | `data/revision/`, `data/parameter_study/` |
| Radius normalization in Appendix A | `data/appendix_figures/` |

[Review notes and plotting instructions](README.txt) identify the data and software associated with each figure. [Environment](environment.json) records the original software versions. Plotting uses additional packages and locally installed Euclid fonts.

## Public copy

The original review archive has SHA-256 `ac1b9a1b64c3efe6f68d4ec60beca34f7b215daf4a73e9e882f7462669791507`. This distribution replaces workstation paths in provenance metadata with relative paths or study identifiers; a historical source identifier need not resolve inside this folder. For MATLAB meshes and example parameters, only the recorded path strings were changed. All numerical arrays were compared with the original archive. The scripts for review and plotting have named steps and portable input paths. The research Python sources are formatted for reading. The 28-feature function was checked against the archived implementation and its saved output. The MATLAB solver entry point documents its inputs and accepted-step outputs. The mesh verifier checks both the public mesh's checksum and the historical checksum used to identify the field simulation. In [SOURCE_INDEX.csv](SOURCE_INDEX.csv), `sha256` identifies the distributed file and `original_sha256` identifies the file in the original review archive. [The file manifest](manifest.json) covers every distributed file apart from the manifest itself.
