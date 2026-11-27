"""Local luminosity correction for a newly inserted column-mass layer."""

from dataclasses import replace
import math

from .spherical_heating import temperature_direction


def correct_inserted_layer(model, index, evaluate, exact_residual, temperature_range, tolerance):
    """Propose a bounded temperature change around one insertion.

    The trial callable returns full-layer sampled luminosity residuals,
    including the recomputed convective flux. A fixed offset calibrates
    those residuals against the initial full-grid transfer.
    """
    if max(map(abs, exact_residual)) <= tolerance:
        return model, 0
    initial = (0.0, 0.0, 0.0)
    controls = range(index - 1, index + 2)
    limit = 0.10
    temperatures = [model.structure[i][3] for i in controls]
    low, high = temperature_range
    bounds = [(max(-limit, low/t - 1), min(limit, high/t - 1)) for t in temperatures]
    trials = {}

    def trial(change):
        if change not in trials:
            rows = list(model.structure)
            for i, delta in zip(controls, change):
                row = rows[i]
                rows[i] = (*row[:3], row[3]*(1 + delta), *row[4:])
            candidate = replace(model, structure=tuple(rows))
            trials[change] = candidate, evaluate(candidate)
        return trials[change]

    sample = trial(initial)[1]
    offset = tuple(full-coarse for full, coarse in zip(exact_residual, sample))

    def errors(change):
        return tuple(a+b for a,b in zip(trial(change)[1], offset))

    def merit(residual):
        return max(map(abs, residual)) + math.sqrt(math.fsum(x*x for x in residual)/len(residual))

    def bounded(values):
        return tuple(max(a, min(b, x)) for x, (a,b) in zip(values, bounds))

    best = initial
    # A warm start probes the convective branch when the remapped layer's
    # local derivative is small because its convective flux is near zero.
    for start in (bounded((0.0, 0.065, 0.035)), initial):
        point = start
        for _ in range(6):
            residual = errors(point)
            if merit(residual) < merit(errors(best)):
                best = point
            if max(map(abs, residual)) <= tolerance:
                break

            columns = []
            for j, (a,b) in enumerate(bounds):
                shift = min(0.001, b-point[j])
                if shift < 0.0005:
                    shift = max(-0.001, a-point[j])
                moved = (*point[:j], point[j]+shift, *point[j+1:])
                columns.append([(x-y)/shift for x,y in zip(errors(moved), residual)])
            direction = temperature_direction(columns, residual)
            for fraction in (1, 0.5, 0.25, 0.125):
                candidate = bounded(x+fraction*d for x,d in zip(point,direction))
                if merit(errors(candidate)) < merit(residual):
                    point = candidate
                    break
            else:
                break

        if merit(errors(point)) < merit(errors(best)):
            best = point
        if max(map(abs, errors(best))) <= tolerance:
            break
    return (model if best == initial else trial(best)[0]), len(trials)
