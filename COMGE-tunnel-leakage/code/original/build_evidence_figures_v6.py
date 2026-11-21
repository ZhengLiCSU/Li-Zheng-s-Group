from pathlib import Path
import sys, json, hashlib, csv, os

for k in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ[k] = "4"
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
import joblib

R = Path(__file__).parents[1]
S = R.parent
O = R / "figures/v6"
O.mkdir(exist_ok=True)
F = (
    S.parent
    / "uq_batches/b96_20261009_rf_reliability_stage01_v3/model/engineering_five_year_20261007/uq_kl_active_20261007"
)
sys.path.insert(0, str(F))
from field_sampler import physical_fields

if os.environ.get("PANEL_AUDIT_DIR"):
    sys.path.insert(0, os.environ["PANEL_AUDIT_DIR"])
from audit_panel_alignment import require_matplotlib_panel_alignment

for p in (S / "paper_manuscript26_20261012/physics/fonts").glob("euclid*.ttf"):
    fm.fontManager.addfont(str(p))
plt.rcParams.update(
    {
        "font.family": "Euclid",
        "font.size": 8.5,
        "axes.labelsize": 8.5,
        "axes.titlesize": 9,
        "legend.fontsize": 7.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "axes.linewidth": 0.6,
        "lines.linewidth": 1.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Euclid",
        "mathtext.it": "Euclid:italic",
        "mathtext.bf": "Euclid:bold",
        "savefig.facecolor": "white",
    }
)
C = ["#1D6483", "#B67530", "#43958B", "#805E99"]
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
outputs = []
bindings = {}


def bind(p):
    bindings[str(p)] = sha(p)
    return p


def label(ax, k, title):
    ax.text(
        0,
        1.07,
        "(" + k + ")",
        transform=ax.transAxes,
        fontweight="bold",
        ha="left",
        va="bottom",
        fontsize=9,
    )
    ax.set_title(title, loc="left", x=0.095, pad=11)
    ax.spines[["top", "right"]].set_visible(False)


def save(fig, name, axes):
    require_matplotlib_panel_alignment(
        fig, axes=list(axes), json_out=O / (name + "_alignment.json"), tolerance_pt=1.5
    )
    fig.savefig(O / (name + ".pdf"))
    fig.savefig(O / (name + ".png"), dpi=600)
    fig.savefig(O / (name + "_preview.png"), dpi=200)
    outputs.extend([name + ".pdf", name + ".png"])
    plt.close(fig)


# Fig. 4: unchanged deterministic histories, visible identities and a peak inset.
M = S / "paper_revision_20261013/mesh_study/postprocess"
cur = pd.read_csv(bind(M / "year5_daily_curves.csv"))
summ = pd.read_csv(bind(M / "mesh_response_summary.csv"))
print("CURVE_COLUMNS", cur.columns.tolist(), flush=True)
days = cur.iloc[:, 0].to_numpy()
vals = [cur.iloc[:, i + 1].to_numpy() for i in range(4)]
fig, aa = plt.subplots(2, 2, figsize=(7.16, 4.7))
fig.subplots_adjust(left=0.11, right=0.98, bottom=0.12, top=0.91, wspace=0.37, hspace=0.65)
ax = aa[0, 0]
styles = ["-", "--", "-.", ":"]
marks = ["o", "s", "^", "D"]
for i, v in enumerate(vals):
    ax.plot(
        days,
        v,
        color=C[i],
        ls=styles[i],
        lw=1.05,
        marker=marks[i],
        mfc="white",
        mew=0.65,
        ms=3,
        markevery=(5 + i * 11, 55),
        label="M" + str(i),
        zorder=5 - i,
    )
ax.set(
    xlabel="Day of year 5",
    ylim=(0, 0.155),
    ylabel=r"Daily leakage ($\mathrm{L\,m^{-1}\,d^{-1}}$)",
    xlim=(0, 365),
)
ax.legend(ncol=4, loc="upper left", columnspacing=0.7, handlelength=1.5, framealpha=1)
label(ax, "a", "Seasonal response")
ins = ax.inset_axes([0.58, 0.13, 0.39, 0.30])
for i, v in enumerate(vals):
    sel = (days >= 180) & (days <= 196)
    ins.plot(days[sel], v[sel], color=C[i], ls=styles[i], lw=0.9)
ins.set(xlim=(181, 195), ylim=(0.1154, 0.1162))
ins.tick_params(labelsize=5.5, pad=1)
ins.set_title("Peak detail", fontsize=6, pad=2, loc="left")
ins.ticklabel_format(useOffset=False)
ins.spines[["top", "right"]].set_visible(False)
ax = aa[0, 1]
for i in range(3):
    ax.plot(
        days,
        100 * (vals[i] - vals[3]) / summ.iloc[3]["year5_annual_L_per_m"],
        color=C[i],
        ls=styles[i],
        label="M" + str(i) + " − M3",
    )
ax.axhline(0, color="#555555", lw=0.5)
ax.set(xlabel="Day of year 5", ylabel="Normalized daily difference (%)", xlim=(0, 365))
ax.legend(frameon=False, fontsize=7)
label(ax, "b", "Signed differences")
ax = aa[1, 0]
v = summ.year5_annual_L_per_m
ax.plot(range(4), v, color="#536778", lw=0.8)
for i, y in enumerate(v):
    ax.scatter(i, y, color=C[i], s=25)
ax.set(
    xticks=range(4),
    xticklabels=["M0", "M1", "M2", "M3"],
    ylabel=r"Year 5 leakage volume ($\mathrm{L\,m^{-1}}$)",
    ylim=(20.875, 20.912),
    xlim=(-0.5, 3.5),
)
label(ax, "c", "Annual leakage")
ax = aa[1, 1]
xs = np.arange(3)
for i, (key, name) in enumerate(
    [
        ("year5_annual_L_per_m_difference_to_M3_percent", "Annual volume"),
        ("five_year_peak24h_L_per_m_day_difference_to_M3_percent", "Five-year rolling peak"),
        ("year5_daily_L1_over_M3_annual_percent", "Daily-history discrepancy"),
    ]
):
    ax.bar(
        xs + (i - 1) * 0.24, abs(summ[key].to_numpy()[:3]), width=0.22, color=C[i], label=name
    )
ax.set(
    xticks=xs,
    xticklabels=["M0", "M1", "M2"],
    ylabel="Absolute difference from M3 (%)",
    ylim=(0, 0.45),
)
ax.legend(frameon=False, fontsize=6.9)
label(ax, "d", "Mesh discrepancies")
save(fig, "fig04_mesh_visible", aa.ravel())

# Fig. 8: actual numerical coordinates; interpolation is exclusively for display.
fld = (
    bind(S / "augmentation128_execution_20261111/fields/nonexistent.npz")
    if False
    else bind(
        S
        / "augmentation128_execution_20261011/adaptive24_design_20261011/fields/adaptive_candidate_000.npz"
    )
)
cfg = json.loads(bind(fld.with_suffix(".json")).read_text("utf8"))["parameters"]
z = np.load(fld)
g = np.load(bind(F / "cache_rf_study_v1/geometry_factors.npz"))
logk, rad, rmean = physical_fields(
    z["soil_gaussian"], z["joint_gaussian"], cfg["soil_cov"], cfg["joint_cov"]
)
xy = np.vstack([g["nodes"], g["quadrature"]])
xx, yy = np.meshgrid(np.linspace(0, 60, 601), np.linspace(0, 40, 401))
cloud = griddata(xy, logk / np.log(10), (xx, yy), method="linear")
cloud[(xx - 30) ** 2 + (yy - 20) ** 2 < 36] = np.nan
cmap = LinearSegmentedColormap.from_list(
    "hydraulic", ["#234D72", "#3D969A", "#E8E8CB", "#C99650", "#994D35"]
)
fig, aa = plt.subplots(1, 2, figsize=(7.16, 2.6))
fig.subplots_adjust(left=0.09, right=0.98, top=0.84, bottom=0.21, wspace=0.40)
ax = aa[0]
im = ax.imshow(
    cloud,
    extent=[0, 60, 40, 0],
    interpolation="none",
    aspect="auto",
    cmap=cmap,
    vmin=-1,
    vmax=1,
)
ax.add_patch(Circle((30, 20), 6, color="#E0E5E8", ec="#36546C", lw=0.5))
ax.add_patch(Circle((30, 20), 5.55, color="white", ec="#36546C", lw=0.5))
ax.set(xlabel="Horizontal coordinate (m)", ylabel="Depth (m)")
label(ax, "a", "Soil hydraulic conductivity")
cbax = ax.inset_axes([0.97, 0.10, 0.025, 0.76])
cb = fig.colorbar(im, cax=cbax)
cb.set_label(r"$\log_{10}(K/\overline{K})$", fontsize=7)
cb.ax.tick_params(labelsize=6)
ax.text(
    0.03,
    0.07,
    "COV = 0.75; SOF = 24 / 6 m",
    transform=ax.transAxes,
    fontsize=6.8,
    bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 2},
)
ax = aa[1]
x = g["joint_coordinates"][:, 0]
ax.plot(x, rad * 1e6, color=C[0], marker="o", ms=2.2, lw=1, label="Realization (64 points)")
ax.axhline(rmean * 1e6, color=C[1], ls="--", label="Ensemble arithmetic mean")
ax.set(xlabel="Distance along J7 pathway (m)", ylabel="Radius scale (µm)", xlim=(0, 0.45))
ax.legend(frameon=False, loc="lower right", fontsize=6.7)
label(ax, "b", "Segment-joint radius scale")
ax.text(0.02, 0.93, "COV = 0.22; SOF = 0.30 m", transform=ax.transAxes, fontsize=7)
save(fig, "fig08_numerical_fields", aa)

# Fig. 12: actual training data and fitted kernel, without neural-network icons.
T = S / "augmentation128_execution_20261011/adaptive24_postprocess_20261011"
train = np.load(bind(T / "training608.npz"))
model = joblib.load(bind(T / "peak608_rbf.joblib"))
X = model["input_scaler"].transform(train["X28"])
gp = model["gp"]
K = gp.kernel_.k1(X)
corr = K / np.sqrt(np.outer(np.diag(K), np.diag(K)))
np.savez_compressed(
    O / "fig12_plot_arrays.npz", standardized_inputs=X, signal_correlations=corr
)
pred = np.load(bind(S / "paper_dense608_prediction_recovery_20261012/mc_results/PC01739.npz"))[
    "predictions"
].reshape(-1)
q = 0.17827224170071532
fig, aa = plt.subplots(2, 2, figsize=(7.16, 4.9))
fig.subplots_adjust(left=0.10, right=0.96, bottom=0.12, top=0.91, wspace=0.41, hspace=0.60)
ax = aa[0, 0]
im = ax.imshow(
    cloud,
    extent=[0, 60, 40, 0],
    interpolation="none",
    aspect="auto",
    cmap=cmap,
    vmin=-1,
    vmax=1,
)
ax.add_patch(Circle((30, 20), 6, fc="white", ec="#36546C", lw=0.6))
ax.set(xlabel="Horizontal coordinate (m)", ylabel="Depth (m)")
label(ax, "a", "Hydraulic input fields")
ins = ax.inset_axes([0.69, 0.61, 0.28, 0.30])
ins.plot(x, rad * 1e6, color=C[0], lw=0.8)
ins.set_title("J7 radius (µm)", fontsize=6, pad=2)
ins.tick_params(labelsize=5.5, pad=1)
ins.spines[["top", "right"]].set_visible(False)
ax = aa[0, 1]
im = ax.imshow(X.T, aspect="auto", vmin=-3, vmax=3, cmap="RdBu_r", interpolation="none")
ax.set(
    xlabel="Training realization",
    ylabel="Standardized input component",
    yticks=[0, 6, 16, 27],
    yticklabels=[1, 7, 17, 28],
)
ax.axhline(6.5, color="white", lw=0.7)
ax.axhline(16.5, color="white", lw=0.7)
label(ax, "b", "28-component representation")
cb = fig.colorbar(im, cax=ax.inset_axes([1.02, 0, 0.025, 1]))
cb.ax.tick_params(labelsize=6)
cb.set_ticks([-3, 0, 3])
ax = aa[1, 0]
im = ax.imshow(
    corr,
    aspect="auto",
    vmin=float(corr.min()),
    vmax=1,
    cmap="viridis",
    interpolation="none",
    origin="lower",
)
ax.set(xlabel="Training realization", ylabel="Training realization")
label(ax, "c", "Fitted Kriging signal correlations")
cb = fig.colorbar(im, cax=ax.inset_axes([1.02, 0, 0.025, 1]))
cb.ax.tick_params(labelsize=6)
cb.set_ticks([round(float(corr.min()), 3), 1])
ax = aa[1, 1]
ax.hist(
    pred,
    bins=55,
    density=True,
    color=C[0],
    alpha=0.5,
    edgecolor="none",
    label="Kriging predictions",
)
ax.axvline(q, color="#B44B46", ls="--", lw=1.3, label="Assessment threshold")
ax.set(
    xlabel=r"Five-year rolling peak ($\mathrm{L\,m^{-1}\,d^{-1}}$)",
    ylabel=r"Probability density ($\mathrm{m\,d\,L^{-1}}$)",
)
label(ax, "d", "Conditional response distribution")
ax.legend(frameon=False, fontsize=6.7)
ax.text(0.97, 0.68, "N = 131,072", transform=ax.transAxes, ha="right", fontsize=7)
for start, end in [
    ((0.465, 0.75), (0.53, 0.75)),
    ((0.725, 0.505), (0.38, 0.505)),
    ((0.49, 0.29), (0.53, 0.29)),
]:
    fig.add_artist(
        FancyArrowPatch(
            start,
            end,
            transform=fig.transFigure,
            arrowstyle="-|>",
            mutation_scale=14,
            lw=1.4,
            color="#263C4A",
        )
    )
save(fig, "fig12_data_workflow", aa.ravel())

# Rebuttal sampling diagnostics from all original case predictions.
v = pd.read_csv(bind(R / "analysis/mc_stability_cases.csv"))
cur = pd.read_csv(bind(R / "analysis/mc_stability_curves.csv"))
pair = pd.read_csv(bind(R / "analysis/mc_stability_paired.csv"))
fig, aa = plt.subplots(2, 2, figsize=(7.16, 4.8))
fig.subplots_adjust(left=0.11, right=0.98, bottom=0.12, top=0.91, wspace=0.40, hspace=0.64)
base = v[v.N == 131072]
ids = [
    base.loc[base.probability.idxmin(), "case_id"],
    base.iloc[(base.probability - 0.05).abs().argmin()].case_id,
    base.loc[base.probability.idxmax(), "case_id"],
]
ax = aa[0, 0]
for i, (cid, lab) in enumerate(
    zip(ids, ["Lowest probability", "Near 5%", "Highest probability"])
):
    t = v[v.case_id == cid]
    ax.plot(t.N / 1024, t.probability * 100, marker="o", ms=3, color=C[i], label=lab)
ax.set(
    xlabel="Total sample size (×1024)",
    ylabel="Exceedance probability (%)",
    xticks=[8, 16, 32, 64, 128],
)
ax.legend(frameon=False, fontsize=6.8, loc="center right", bbox_to_anchor=(1, 0.70))
label(ax, "a", "Sample-size sensitivity")
ax = aa[0, 1]
a = v[v.N == 65536].set_index("case_id")
b = v[v.N == 131072].set_index("case_id")
dp = (b.probability - a.probability) * 100
ax.hist(dp, bins=30, color=C[0], alpha=0.75, edgecolor="white", lw=0.2)
ax.axvline(0, color="#555555", lw=0.6, ymax=0.8)
ax.set(xlabel="Probability change (percentage points)", ylabel="Number of parameter settings")
label(ax, "b", "65,536 to 131,072 samples")
ax.text(
    0.97, 0.93, "All 4,181 settings", transform=ax.transAxes, va="top", ha="right", fontsize=7
)
ax = aa[1, 0]
for i, (cid, lab) in enumerate(
    zip(ids, ["Lowest probability", "Near 5%", "Highest probability"])
):
    t = v[v.case_id == cid]
    ax.plot(t.N / 1024, t.p95, marker="o", ms=3, color=C[i], label=lab)
ax.set(
    xlabel="Total sample size (×1024)",
    ylabel=r"95th percentile ($\mathrm{L\,m^{-1}\,d^{-1}}$)",
    xticks=[8, 16, 32, 64, 128],
)
label(ax, "c", "Response-quantile sensitivity")
ax = aa[1, 1]
selected = pair[
    (pair.N == 131072)
    & pair.group_id.isin(
        [
            "soil_horizontal_exponential_cv0.5",
            "soil_vertical_exponential_cv0.5",
            "soil_proportional_exponential_cv0.5",
        ]
    )
]
labs = []
for i, row in enumerate(selected.itertuples()):
    ax.errorbar(
        row.delta_p * 100, i, xerr=row.delta_p_SE * 100, fmt="o", ms=4, capsize=3, color=C[i]
    )
    labs.append(row.axis.capitalize())
ax.axvline(0, color="#777777", lw=0.7, ymax=0.8)
ax.set(
    yticks=range(len(labs)),
    yticklabels=labs,
    xlabel="Paired probability difference (pp)",
    ylim=(-0.6, len(labs) - 0.4),
)
ax.text(
    0.03,
    0.96,
    "Error bars: ±1 replicate standard error",
    transform=ax.transAxes,
    fontsize=6.6,
    va="top",
)
label(ax, "d", "Small soil-correlation effects")
save(fig, "response_mc_stability", aa.ravel())
with (O / "figure_provenance.json").open("x", encoding="utf8") as f:
    json.dump(
        {
            "passed": True,
            "source_sha256": sha(Path(__file__)),
            "bindings": bindings,
            "outputs": {n: sha(O / n) for n in outputs},
            "physical_solver_calls": 0,
            "new_fits": 0,
            "AI_generated_content": False,
            "font": "Euclid",
        },
        f,
        indent=2,
        ensure_ascii=False,
    )
print("FIGURES_COMPLETED", outputs, flush=True)
