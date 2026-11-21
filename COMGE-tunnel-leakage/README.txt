REVIEW REPRODUCIBILITY ARCHIVE
Tunnel leakage reliability under freeze-thaw cycling

This archive accompanies the revised manuscript and its response to Reviewer 2.
It provides exact feature definitions, a stored full field, the frozen Kriging
predictor, training and held-out data, split records, selection rules, numerical
controls, source files and the saved statistics cited in the response.

Start here
1. Extract the complete archive into a new directory.
2. Use Python with NumPy and SciPy. environment.json lists the versions used in
   the successful review check. The portable predictor does not load joblib or
   require scikit-learn. The original joblib model is included for provenance.
3. From the extracted directory run
       python verify_review.py --output-dir review_check
   Choose a new output directory if that name already exists.
4. Read review_check/verification.json. The program also exports all 80 replayed
   validation predictions and the exact histogram-bin data for Fig.12(d).

The check reconstructs the 28 features from the stored full numerical field and
the original feature function. It then evaluates the frozen predictor through
its explicit fitted arrays and compares all 80 saved held-out predictions.
It verifies case separation and both CV split schemes, all 20,905 sample-size
records, 270 paired contrasts and the reported paired mesh differences.
It does not launch MATLAB, fit a model, sample new fields or modify the source
data. These operations are a reproducibility check of completed work.

Navigation by manuscript and reviewer concern
Section 2.2 / Comment 4
  solver/source/ contains the unmodified custom MATLAB files.
  solver/numerical_controls.json records executed solver, integration, initial
  and mechanical settings. solver/example/executed_parameters.mat contains
  the complete executed parameter structure. solver/meshes/ gives M0 and M1.
  solver/source_context/run_matched_case.m documents the original invocation.
  Its original once/observer controls are historical execution safeguards, not
  a portable restart command. The review entry point does not execute it.

Section 2.3 / Comment 6
  data/revision/matched_physics_results.json gives all six completed responses,
  peak definitions and paired comparisons. The verifier independently integrates
  all six interval-rate histories and checks the peaks and fifth-year volumes. The RF*_response.npz files contain
  their saved histories. This comparison changes the soil far field. The local
  lining/interface mesh and 64-point segment-joint profile are fixed. It does
  not test heterogeneous local/J7 refinement or full probability convergence.

Sections 3.1 and 4.2 / Comments 7-12
  code/original/field_sampler.py, hybrid_kl.py and uq_workflow.py retain the exact
  generation and feature definitions. The original scalar feature routine is
  called directly by verify_review.py.
  data/example/geometry.npz contains all 48,479 soil coordinates, witness indices
  and the 64 segment-joint coordinates. soil_feature_sampling.csv provides the
  exact 169 coordinates, regional memberships and equal weights.
  joint_feature_sampling.csv defines four consecutive blocks of 16 points.
  data/example/adaptive_candidate_000.npz contains the physical field shown in Fig. 8.
  Fig. 12(a) uses the first Monte Carlo realization in data/figure12/PC01739.npz.
  This field identifier is an internal file key, not a new experimental result.
  FEATURE_DEFINITIONS.txt gives the ordered 28-component representation.

Section 4.2 / Comment 12
  data/model/training608.npz contains the ordered training inputs and responses.
  case_roles.json distinguishes the original 512, prespecified 72, adaptive 24
  and the 80 held-out cases. training608_evaluation.json gives every fold's
  exact training/test IDs and parameter groups, scaling and fitted kernel.
  data/selection/ contains the recorded candidate scores, selected inputs,
  family diagnostics and original selection code.
  data/model/frozen_predictor_arrays.npz is a lossless numerical export of the
  original fitted model, with input scaling, training coordinates, kernel
  parameters, solved coefficients and response scaling.
  data/validation/ contains all held-out data and threshold diagnostics.
  These response checks include two exceedances overall and one exceedance in
  the eleven-case subset within 20% of the threshold. They do not calibrate
  physical probability error over the entire parameter domain.

Fig.12(d)
  The original NumPy generator used SeedSequence([20261012,60821,rep,stream])
  with rep = 0 through 7, stream 0 for soil and stream 1 for the segment joint.
  data/revision/gaussian_stream_receipt.json records both Gaussian array hashes.
  data/figure12/PC01739.npz contains the 131,072 actual predictions used.
  The adjacent receipt records every parameter at full precision.
  Soil COV 0.25, segment-joint COV 0.25, soil SOFs 25.722563100871 and
  6.430640775218 m, segment-joint SOF 0.30 m. Soil autocorrelation is squared
  exponential and segment-joint autocorrelation is exponential. F4 normalization
  and q*=0.17827224170071532 L m^-1 d^-1 apply.
  Panel (a) shows the first realization at this setting. Panel (c) predicts
  its peak, and panel (d) uses all 131,072 predictions at the same setting.

Sections 4.3-4.6 / Comment 13
  data/revision/mc_stability_* contains all nested sample-size statistics and
  eight-replication arrays. Probabilities are fractions unless explicitly
  labelled percentage points. Replicate contrasts preserve common random
  numbers. Sample-size relative mean and quantile changes use the final
  131,072-sample estimate as denominator.
  data/parameter_study/ supplies the full processed tables for all 4,181
  parameter settings, 54 curves with 26 evaluated points each, five 26x26
  surfaces, the 41-threshold sweep, paired contrasts and response perturbations.
  Full raw predictions for every setting and every physical transient history
  are not duplicated here. The selected raw distribution, processed tables
  and replay records make the checks above portable; the complete source
  archive remains available from the corresponding author.

Provenance and interpretation
manifest.json binds every archive file except itself, including SOURCE_INDEX.csv.
SOURCE_INDEX.csv records the purpose and original SHA256 of copied material.
The submitted manuscript and response remain the authority for physical scope.
Two new five-year physical runs are supplied for the local mesh comparison in
response Fig. R2.2. The main reliability results use the original fitted model.
F2 probability shifts remain fixed-surrogate normalization sensitivities; the
three F2 physical checks are below threshold and do not validate that tail.


Revised Fig. 12 and response Fig. R2.5
From the extracted archive, run python code/plot_figure12.py --output figure12
with NumPy, SciPy, Matplotlib and Euclid regular/bold fonts. If Euclid is
not installed, pass --font-dir pointing to the directory containing the fonts.
Panel (a) plots the first realization stored in data/figure12/PC01739.npz
at the 169 soil coordinates and 64 J7 coordinates used by the predictor.
Panel (b) pairs each training input with its THM peak. Panel (c) reconstructs
the new input and its kernel correlation to every training realization,
then reproduces the saved peak through the fitted coefficients. Panel (d)
uses the original predictions and 55 histogram bins. The count is taken
directly from the predictions above the threshold. This script supersedes
the Fig. 12 section of code/original/build_evidence_figures_v6.py.

Figures 8 and A.1-A.2 (also R2.4 and R3.1)
Run python code/plot_field_and_appendix.py --output field_and_appendix from this extracted folder.
The figure script requires NumPy, SciPy, Matplotlib, Pillow, PyMuPDF and Euclid fonts.
Pass --font-dir for locally installed Euclid TrueType files if needed.
Figure 8 uses the saved soil and segment-joint realization in data/example. Appendix vector sources and numeric comparison tables are in data/appendix_figures.

Local refinement near the J7 threshold (response Fig. R2.2)
The two cases were selected by the nearest physical peak on either side of the
threshold in the original 688 simulations. The refined mesh retains the far
field and reduces the J7 band spacing from 10 to 5 mm. The two lines below
recalculate the peaks from their accepted-step rates and reproduce the figure.

    python code/verify_threshold_mesh.py
    python code/plot_threshold_mesh.py --data data/revision/near_threshold_mesh --output Fig_R2_2.png

The rates in rolling_rate_histories.npz are per metre of tunnel axis in litres
per metre per day. Each row stores the end time (days) and the rate for the
preceding accepted interval. selection.csv gives the original physical peaks.
peaks.csv contains the peak, the start of the maximizing 24-hour window and
the percentage margin 100 times (peak/threshold minus one). The mesh view
contains the actual element vertices in millimetres about the outlet. The
figure script uses Euclid fonts, with --fonts for an alternative font directory.
The MATLAB files contain the original and locally refined mesh geometries and
the hydraulic field inputs for both realizations on both meshes. Soil
coordinates are in metres; log multipliers are dimensionless. The segment-joint
radius profiles are in metres. The verification uses NumPy and SciPy; the figure also
uses Matplotlib and Pillow.
