# Tunnel leakage reliability under freeze–thaw cycling

Companion data and code for the *Computers and Geotechnics* manuscript on the reliability of leakage through a prescribed tunnel segment joint. The code reconstructs the random field inputs, replays the fitted Kriging predictions, and checks the Monte Carlo and mesh-sensitivity results discussed in the response to Reviewer 2.

[Code guide](CODE_GUIDE.md) · [Kriging verification](verify_review.py) · [Hydraulic feature construction](code/original/uq_workflow.py)

## Reproduce the results

Download the repository or clone it. From **this directory**, use Python to reconstruct the complete study directory at a fresh location. Reassembly uses the Python standard library. It checks every reconstructed file against the scientific package's SHA-256 manifest.

```sh
python -B prepare_review.py --output-dir ../tunnel_leakage_study
cd ../tunnel_leakage_study
python -B verify_review.py --output-dir ../tunnel_leakage_review_results
python -B code/verify_threshold_mesh.py
```

The paths in the code guide refer to the reconstructed study directory. The two verification commands require NumPy and SciPy. The first writes the reconstruction results and checks the predictions on 80 held-out cases. The second re-integrates discharge histories for the mesh comparisons on either side of the leakage threshold.

## Read the data directly

[Study overview and figure guide](delivery_parts/README.study.md) · [Detailed review notes](README.txt) · [Software versions](environment.json) · [Input definitions](FEATURE_DEFINITIONS.txt).

Large CSV tables are saved as numbered CSV files under [delivery_parts](delivery_parts/); each has the full original header and a consecutive, non-overlapping set of data rows. The [partition index](delivery_parts/index.json) lists the source filename, ordered pieces, row counts for CSVs, sizes and checksums. Model and NumPy binary pieces are joined in that order by the preparation script. All other scientific files are provided directly. Reassembly recovers the original full CSVs, model and data files byte for byte.
