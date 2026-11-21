# Code guide

## From the hydraulic fields to a leakage prediction

[verify_review.py](verify_review.py) starts with a saved realization of the soil conductivity and segment-joint radius fields. It reads the geometry, Gaussian values and covariance parameters from [data/example](data/example), then calls `physical_features` in [uq_workflow.py](code/original/uq_workflow.py). This produces the 28 inputs used by the surrogate. [FEATURE_DEFINITIONS.txt](FEATURE_DEFINITIONS.txt) lists them in model order.

To reproduce a prediction, `predict` in the verifier standardizes the inputs using the training means and scales. It evaluates the squared-exponential kernel against the 608 training inputs, multiplies these covariances by the fitted coefficients, then reverses the scaling and log transformation of the leakage response. All fitted quantities are stored in [frozen_predictor_arrays.npz](data/model/frozen_predictor_arrays.npz). The verifier compares the predictions with the saved predictions for 80 held-out cases. [plot_figure12.py](code/plot_figure12.py) illustrates this calculation for an individual draw, then plots the distribution for all draws at that setting.

## Coupled tunnel model

The five-year simulations used [TH_2D_segment_v16.m](solver/source/TH_2D_segment_v16.m). Its input `opt` specifies the mesh, properties, initial conditions and numerical controls. The returned structure `out` contains the accepted time-step history and the leakage rates. These functions locate the main calculations within the implementation.

| Calculation | MATLAB source |
| --- | --- |
| Soil integration and water flow | `build_soil_scv_quadrature.m`, `build_soil_cvfe_flux.m` and `resolve_soil_conductivity_field.m` |
| Coupled soil response and lining mesh | `thm_bulk_constitutive.m` and `build_segment_mechanical_mesh.m` |
| Segment-joint conductance, channel storage and drainage | `joint_freeze_conductance.m`, `joint_channel_inventory.m` and `evaluate_joint_drain.m` |

These functions are in [solver/source](solver/source). The matched-mesh example supplies the [executed parameters and input fields](solver/example), the [meshes](solver/meshes) and the [original MATLAB invocation](solver/source_context/run_matched_case.m). That invocation expects the original execution folders, case plan and process observer. The review scripts use the stored discharge histories to recalculate the reported peaks.

## Running the checks

From the directory containing this guide, run the following commands with Python, NumPy and SciPy. Choose a fresh output directory for the first command.

```sh
python -B verify_review.py --output-dir ../tunnel_review_results
python -B code/verify_threshold_mesh.py
```

The first command writes `verification.json`, `validation_replayed.csv` and `figure12d_histogram.csv`. The second recomputes the rolling 24-hour peak from accepted flow intervals on the original and refined meshes for the physical cases immediately below and above the threshold. [The review notes](README.txt) give the plotting commands for Figs. 8, 12, A.1 and A.2 and the mesh comparison in the response.

## Training and adaptive selection

[fit_peak608_once.py](code/original/fit_peak608_once.py) records how the original calculation applied grouped cross-validation and trained the final Kriging model. [The selection records](data/selection) identify the candidates, their scores and the 24 cases selected for supplementary training. The original campaign scripts refer to intermediate model fits, campaign manifests and helper modules stored with that campaign. The portable verification uses the recorded splits, selection results and fitted model arrays.
