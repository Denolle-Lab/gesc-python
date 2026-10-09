"""Consistent scientific figures and machine-readable notebook checks."""

import html
import json

import matplotlib.pyplot as plt
import numpy as np
from IPython.display import HTML, display

COLORS = ("#0f766e", "#d97706", "#7c3aed", "#b91c1c")


def style():
    plt.rcParams.update(
        {
            "figure.dpi": 130,
            "savefig.dpi": 160,
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.prop_cycle": plt.cycler(color=COLORS),
            "axes.grid": True,
            "grid.alpha": 0.18,
            "figure.constrained_layout.use": True,
        }
    )


def geological_labels(operator):
    return [
        f"{p.top_m:g}–{p.bottom_m:g} m" if p.bottom_m is not None else f"{p.top_m:g} m–∞"
        for p in operator.parameters
    ]


def model_table(model):
    rows = [
        (f"{r.top_m:g}–{r.bottom_m:g}", r.top.vp_m_s, r.top.vs_m_s, r.top.density_kg_m3)
        for r in model.regions
    ] + [
        (
            f"{model.halfspace_top_m:g}–∞",
            model.halfspace.vp_m_s,
            model.halfspace.vs_m_s,
            model.halfspace.density_kg_m3,
        )
    ]
    body = "".join(f"<tr>{''.join(f'<td>{value}</td>' for value in row)}</tr>" for row in rows)
    display(
        HTML(
            "<table><thead><tr><th>Depth (m)</th><th>Vp (m/s)</th>"
            "<th>Vs (m/s)</th><th>Density (kg/m³)</th></tr></thead><tbody>"
            + body
            + "</tbody></table>"
        )
    )


def sensitivity_figure(modes, operator):
    frequency = np.array(operator.observation.frequencies_hz)
    labels = geological_labels(operator)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7))
    axes[0].plot(frequency, [m.phase_velocity_m_s for m in modes], "o-")
    axes[0].set(
        xlabel="Frequency (Hz)", ylabel="Phase velocity (m/s)", title="Fundamental Rayleigh"
    )
    for j, label in enumerate(labels):
        axes[1].plot(frequency, np.array(operator.g)[:, j], "o-", label=label)
    axes[1].set(
        xlabel="Frequency (Hz)",
        ylabel="∂ ln c / ∂ ln Vs (region)",
        title="Integrated layer sensitivity",
    )
    axes[1].legend(fontsize=8)
    for mode in (modes[0], modes[-1]):
        depth = np.array(mode.depth_m)
        amplitude = np.hypot(mode.horizontal, mode.vertical)
        axes[2].plot(amplitude, depth, label=f"{mode.frequency_hz:g} Hz")
    limit = max(operator.parameters[-1].top_m * 1.5, 80)
    axes[2].set(
        xlabel="Relative displacement",
        ylabel="Depth (m)",
        ylim=(limit, 0),
        title="Eigenfields: not sensitivity kernels",
    )
    axes[2].legend(fontsize=8)
    fig.suptitle("Synthetic structure • frequency samples overlapping depths", fontsize=12)
    plt.show()


def report_figure(operator, truth, report, title):
    mean = np.array(report.mean_log_vs)
    sd = np.sqrt(np.diag(report.covariance_log_vs))
    labels = geological_labels(operator)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    pos = np.arange(len(mean))
    axes[0].errorbar(mean * 100, pos, xerr=1.96 * sd * 100, fmt="o", label="Posterior 95%")
    axes[0].plot(np.array(truth) * 100, pos, "x", color=COLORS[1], label="Injected truth")
    axes[0].set(yticks=pos, yticklabels=labels, xlabel="Δ ln Vs (%)", title="Depth coefficients")
    axes[0].invert_yaxis()
    axes[0].legend(fontsize=8)
    averaging = np.array(report.averaging_matrix)
    image = axes[1].imshow(averaging, cmap="RdBu_r", vmin=-1, vmax=1)
    axes[1].set(
        xticks=pos,
        yticks=pos,
        xticklabels=labels,
        yticklabels=labels,
        xlabel="True coefficient",
        ylabel="Estimated coefficient",
        title="Averaging matrix",
    )
    axes[1].tick_params(axis="x", rotation=25)
    fig.colorbar(image, ax=axes[1], shrink=0.7)
    correlation = np.array(report.covariance_log_vs) / np.outer(sd, sd)
    image = axes[2].imshow(correlation, cmap="RdBu_r", vmin=-1, vmax=1)
    axes[2].set(
        xticks=pos,
        yticks=pos,
        xticklabels=labels,
        yticklabels=labels,
        title="Posterior correlation",
    )
    axes[2].tick_params(axis="x", rotation=25)
    fig.colorbar(image, ax=axes[2], shrink=0.7)
    fig.suptitle(title)
    plt.show()


def emit_metrics(metrics):
    display({"application/json": metrics, "text/plain": json.dumps(metrics, indent=2)}, raw=True)
    rows = "".join(
        f"<tr><td>{html.escape(key.replace('_', ' '))}</td><td>{html.escape(str(value))}</td></tr>"
        for key, value in metrics.items()
    )
    display(HTML("<table><caption>Executed checks and diagnostics</caption>" + rows + "</table>"))
