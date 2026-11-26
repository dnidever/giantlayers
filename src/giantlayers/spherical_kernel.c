/* Optional compiled wavelength kernel for spherical coherent scattering. */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#if defined(_WIN32) && defined(GIANTLAYERS_EXTENSION)
#define KERNEL_EXPORT __declspec(dllexport)
#else
#define KERNEL_EXPORT
#endif

static void moments(int n, int rays, const int *offset, const int *core,
                    const double *transmit, const double *start, const double *end,
                    const double *bottom, const double *jw, const double *fw,
                    const double *source, double *mean, double *flux, double *incoming) {
    for (int i = 0; i < n; ++i) {
        mean[i] = 0.0;
        if (flux) flux[i] = 0.0;
    }
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

static int direct_scattering(int n, int rays, const int *offset, const int *core,
                             const double *transmit, const double *start, const double *end,
                             const double *bottom, const double *jw, const double *fw,
                             const double *thermal, const double *fraction, double *source,
                             double *incoming) {
    double *matrix = malloc(((size_t)n + 4) * n * sizeof(*matrix));
    if (!matrix) return -2;
    double *rhs = matrix + (size_t)n*n;
    double *background = rhs+n, *basis = background+n, *response = basis+n;
    for (int i = 0; i < n; ++i) basis[i] = 0.0;
    moments(n, rays, offset, core, transmit, start, end, bottom, jw, fw, basis, background, NULL, incoming);
    for (int column = 0; column < n; ++column) {
        basis[column] = 1.0;
        moments(n, rays, offset, core, transmit, start, end, bottom, jw, fw, basis, response, NULL, incoming);
        for (int row = 0; row < n; ++row)
            matrix[(size_t)row*n+column] = (row == column) - fraction[row]*(response[row]-background[row]);
        basis[column] = 0.0;
    }
    for (int i = 0; i < n; ++i)
        rhs[i] = (1.0-fraction[i])*thermal[i] + fraction[i]*background[i];

    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int row = i+1; row < n; ++row)
            if (fabs(matrix[(size_t)row*n+i]) > fabs(matrix[(size_t)pivot*n+i])) pivot = row;
        if (fabs(matrix[(size_t)pivot*n+i]) < 1e-15) {
            free(matrix);
            return -1;
        }
        if (pivot != i) {
            double swap = rhs[i];
            rhs[i] = rhs[pivot];
            rhs[pivot] = swap;
            for (int j = i; j < n; ++j) {
                swap = matrix[(size_t)i*n+j];
                matrix[(size_t)i*n+j] = matrix[(size_t)pivot*n+j];
                matrix[(size_t)pivot*n+j] = swap;
            }
        }
        for (int row = i+1; row < n; ++row) {
            double factor = matrix[(size_t)row*n+i]/matrix[(size_t)i*n+i];
            rhs[row] -= factor*rhs[i];
            for (int j = i+1; j < n; ++j)
                matrix[(size_t)row*n+j] -= factor*matrix[(size_t)i*n+j];
        }
    }
    for (int i = n-1; i >= 0; --i) {
        double value = rhs[i];
        for (int j = i+1; j < n; ++j) value -= matrix[(size_t)i*n+j]*source[j];
        source[i] = value/matrix[(size_t)i*n+i];
    }
    free(matrix);
    return 0;
}

KERNEL_EXPORT int spherical_solve(int n, int rays, const int *offset, const int *core,
                    const double *scale, const double *mu_bottom,
                    const double *jw, const double *fw, const double *radial,
                    const double *thermal, const double *absorption, const double *scattering,
                    double tolerance, int max_iterations, double *mean, double *flux, double *heating) {
    int entries = offset[rays];
    double *buffers = malloc((3*(size_t)entries + (size_t)rays + 5*(size_t)n) * sizeof(*buffers));
    if (!buffers) return -2;
    double *transmit = buffers, *start = transmit+entries, *end = start+entries, *bottom = end+entries;
    double *diagonal = bottom+rays, *source = diagonal+n, *fraction = source+n;
    double *updated = fraction+n, *incoming = updated+n;
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
        if (core[ray] && radial[n-2] > 0) {
            double x = fmax(0.0, fmin(1.0, (distance-5.0)/10.0));
            double increment = mu_bottom[ray]*(thermal[n-1]-thermal[n-2])/radial[n-2];
            bottom[ray] = fmax(0.0, bottom[ray] + x*x*(3.0-2.0*x)*increment);
        }
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
            moments(n, rays, offset, core, transmit, start, end, bottom, jw, fw, source, mean, NULL, incoming);
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
        if (iterations > max_iterations) {
            int status = direct_scattering(n, rays, offset, core, transmit, start, end, bottom,
                                           jw, fw, thermal, fraction, source, incoming);
            if (status < 0) {
                free(buffers);
                return status;
            }
            iterations = max_iterations;
        }
    }
    moments(n, rays, offset, core, transmit, start, end, bottom, jw, fw, source, mean, flux, incoming);
    for (int i = 0; i < n; ++i)
        heating[i] = 12.56637061435917295385*absorption[i]*(mean[i]-thermal[i]);
    free(buffers);
    return iterations;
}
