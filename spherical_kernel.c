/* Optional compiled wavelength kernel for spherical coherent scattering. */
#include <math.h>
#include <stddef.h>

static void moments(int n, int rays, const int *offset, const int *core,
                    const double *transmit, const double *start, const double *end,
                    const double *bottom, const double *jw, const double *fw,
                    const double *source, double *mean, double *flux) {
    for (int i = 0; i < n; ++i) {
        mean[i] = 0.0;
        if (flux) flux[i] = 0.0;
    }
    double incoming[n];
    for (int ray = 0; ray < rays; ++ray) {
        int base = offset[ray], last = offset[ray + 1] - base - 1;
        incoming[0] = 0.0;
        for (int i = 0; i < last; ++i) {
            int k = base + i;
            incoming[i+1] = transmit[k]*incoming[i] + start[k]*source[i] + end[k]*source[i+1];
        }
        double outgoing = core[ray] ? bottom[ray] : incoming[last];
        for (int i = last; i >= 0; --i) {
            int k = base + i;
            if (i < last)
                outgoing = transmit[k]*outgoing + start[k]*source[i+1] + end[k]*source[i];
            mean[i] += jw[k]*(outgoing+incoming[i]);
            if (flux) flux[i] += fw[k]*(outgoing-incoming[i]);
        }
    }
}

int spherical_solve(int n, int rays, const int *offset, const int *core,
                    const double *scale, const double *mu_bottom,
                    const double *jw, const double *fw, const double *radial,
                    const double *thermal, const double *absorption, const double *scattering,
                    double tolerance, int max_iterations, double *mean, double *flux, double *heating) {
    int entries = offset[rays];
    double transmit[entries], start[entries], end[entries], bottom[rays];
    double diagonal[n], source[n], fraction[n], updated[n];
    double peak = 0.0;
    int has_scattering = 0;
    for (int i = 0; i < n; ++i) {
        double opacity = absorption[i] + scattering[i];
        fraction[i] = opacity ? scattering[i]/opacity : 0.0;
        if (fraction[i]) has_scattering = 1;
        diagonal[i] = 0.0;
        source[i] = thermal[i];
        peak = fmax(peak, thermal[i]);
    }
    for (int ray = 0; ray < rays; ++ray) {
        int base = offset[ray], last = offset[ray+1] - base - 1;
        double distance = 0.0;
        for (int i = 0; i < last; ++i) {
            int k = base+i;
            double step = radial[i]*scale[k];
            distance += step;
            transmit[k] = exp(-step);
            double absorbed = -expm1(-step);
            end[k] = step < 1e-4 ? step*(0.5-step/6.0+step*step/24.0) : 1-absorbed/step;
            start[k] = absorbed-end[k];
        }
        bottom[ray] = thermal[n-1];
        if (core[ray] && distance >= 10 && radial[n-2] > 0)
            bottom[ray] += mu_bottom[ray]*(thermal[n-1]-thermal[n-2])/radial[n-2];
        double tangent_transmission = 1.0;
        for (int i = last; i >= 0; --i) {
            int k = base + i;
            double inbound = i ? end[k-1] : 0.0;
            double outbound = i < last ? end[k] : 0.0;
            if (!core[ray]) {
                if (i == last) {
                    outbound = inbound;
                } else {
                    outbound += start[k]*transmit[k]*tangent_transmission;
                    tangent_transmission *= transmit[k]*transmit[k];
                    outbound += inbound*tangent_transmission;
                }
            }
            diagonal[i] += jw[k]*(inbound+outbound);
        }
    }

    int iterations = 0;
    if (has_scattering && peak > 0) {
        for (iterations = 1; iterations <= max_iterations; ++iterations) {
            moments(n, rays, offset, core, transmit, start, end, bottom, jw, fw, source, mean, NULL);
            double error = 0.0;
            for (int i = 0; i < n; ++i) {
                double w = fraction[i];
                double residual = (1-w)*thermal[i] + w*mean[i] - source[i];
                updated[i] = source[i] + residual/(1-w*diagonal[i]);
                error = fmax(error, fabs(updated[i]-source[i])/fmax(fabs(updated[i]), peak*1e-12));
            }
            for (int i = 0; i < n; ++i) source[i] = updated[i];
            if (error <= tolerance) break;
        }
        if (iterations > max_iterations) return -1;
    }
    moments(n, rays, offset, core, transmit, start, end, bottom, jw, fw, source, mean, flux);
    for (int i = 0; i < n; ++i)
        heating[i] = 12.56637061435917295385*absorption[i]*(mean[i]-thermal[i]);
    return iterations;
}
