"""Plot the structure and luminosity diagnostics of a spherical atmosphere."""

from math import log10
from pathlib import Path
from tempfile import NamedTemporaryFile

PLOT_FORMATS = {".png", ".pdf", ".svg"}


def write_profile_plot(path, model, depth, budget, *, teff, logg, status, maximum, tolerance):
    """Save the final atmosphere profiles against Rosseland optical depth."""
    from matplotlib.figure import Figure

    log_depth = [log10(value) if value > 0 else float("nan") for value in depth]
    temperature = [row[3] for row in model.structure]
    gas_pressure = [row[5] for row in model.structure]
    electron_pressure = [row[4] for row in model.structure]
    density = [row[2] for row in model.thermodynamics]
    residual = [100 * row[-1] for row in budget]

    fig = Figure(figsize=(9, 6.5), layout="constrained")
    (temp, pressure), (rho, error) = fig.subplots(2, 2, sharex=True)
    temp.plot(log_depth, temperature, color="#1d577e", marker=".", linewidth=1.5)
    temp.set_ylabel("Temperature (K)")
    pressure.semilogy(log_depth, gas_pressure, color="#1d577e", marker=".", label="Gas")
    pressure.semilogy(log_depth, electron_pressure, color="#a44d20", linestyle="--", label="Electron")
    pressure.set_ylabel(r"Pressure (dyn cm$^{-2}$)")
    pressure.legend(frameon=False)
    rho.semilogy(log_depth, density, color="#1d577e", marker=".")
    rho.set_ylabel(r"Mass density (g cm$^{-3}$)")
    error.axhspan(-100 * tolerance, 100 * tolerance, color="#cbe5db", alpha=0.7)
    error.axhline(0, color="#777777", linewidth=0.8)
    error.plot(log_depth, residual, color="#1d577e", marker=".")
    error.set_ylabel("Luminosity residual (%)")

    for ax in (temp, pressure, rho, error):
        if depth[0] <= 1 <= depth[-1]:
            ax.axvline(0, color="#858585", linestyle=":", linewidth=1)
        ax.grid(alpha=0.2)
    for ax in (rho, error):
        ax.set_xlabel(r"$\log_{10}\,\tau_\mathrm{Ross}$")
    fig.suptitle(f"GiantLayers  |  {teff:g} K, log g = {logg:g}\n"
                 f"{status.capitalize()} · peak luminosity residual {maximum:.2%}")

    with NamedTemporaryFile(dir=path.parent, suffix=path.suffix, delete=False) as temporary:
        target = Path(temporary.name)
    try:
        fig.savefig(target, format=path.suffix.lower().lstrip("."), dpi=160)
        target.replace(path)
    finally:
        target.unlink(missing_ok=True)
