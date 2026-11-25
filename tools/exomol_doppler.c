#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Offline ExoMol line integration. All spectral coordinates are in cm^-1. */
static int containing_bin(const double *edges, int bins, double value)
{
    int low = 0, high = bins;
    while (high - low > 1) {
        int mid = low + (high - low)/2;
        if (edges[mid] > value) high = mid;
        else low = mid;
    }
    return low;
}
int main(int argc, char **argv)
{
    if (argc != 2) return 1;
    FILE *config = fopen(argv[1], "r");
    if (!config) return 1;
    int bins, nt, states;
    double mass, vmic;
    if (fscanf(config, "%d%d%d%lf%lf", &bins, &nt, &states, &mass, &vmic) != 5 ||
        bins < 1 || nt < 1 || states < 1 || !(mass > 0) || !(vmic >= 0)) return 1;
    double *edges = calloc((size_t)bins+1, sizeof(double));
    double *speed = calloc(nt, sizeof(double));
    double *energy = calloc(states+1, sizeof(double));
    double *weight = calloc(states+1, sizeof(double));
    double *t = calloc(nt, sizeof(double));
    double *q = calloc(nt, sizeof(double));
    double *spectrum = calloc((size_t)nt*bins, sizeof(double));
    double *boltzmann = calloc((size_t)nt*(states+1), sizeof(double));
    if (!edges || !speed || !energy || !weight || !t || !q || !spectrum || !boltzmann) return 1;
    const double c = 2.99792458e10, c2 = 1.438776877, kb = 1.380649e-16;
    const double amu = 1.66053906660e-24, pi = 3.14159265358979323846;
    for (int i=0; i<=bins; ++i)
        if (fscanf(config, "%lf", &edges[i]) != 1 || !(edges[i] > (i ? edges[i-1] : 0))) return 1;
    double max_speed = 0;
    for (int j=0; j<nt; ++j) {
        if (fscanf(config, "%lf%lf", &t[j], &q[j]) != 2 || !(t[j] > 0) || !(q[j] > 0)) return 1;
        speed[j] = sqrt(2*kb*t[j]/(mass*amu) + vmic*vmic)/c;
        max_speed = fmax(max_speed, speed[j]);
    }
    int id;
    double e, g;
    while (fscanf(config, "%d%lf%lf", &id, &e, &g) == 3) {
        if (id < 1 || id > states || weight[id] != 0 || !isfinite(e) || !(g > 0)) return 1;
        energy[id] = e;
        weight[id] = g;
        for (int j=0; j<nt; ++j)
            boltzmann[(size_t)id*nt+j] = exp(-c2*e/t[j])/q[j];
    }
    fclose(config);
    char line[512];
    unsigned long count=0, used=0;
    while (fgets(line, sizeof(line), stdin)) {
        char *end;
        long upper = strtol(line, &end, 10);
        long lower = strtol(end, &end, 10);
        double A = strtod(end, &end);
        if (upper < 1 || lower < 1 || upper > states || lower > states ||
            weight[upper] == 0 || weight[lower] == 0 || !(A >= 0)) return 2;
        ++count;
        double wn = energy[upper] - energy[lower];
        if (wn <= 0 || wn*(1+5*max_speed) < edges[0] ||
            wn*(1-5*max_speed) > edges[bins] || A == 0) continue;
        ++used;
        double factor = A*weight[upper]/(8*pi*c*wn*wn);
        int center = containing_bin(edges, bins, wn);
        for (int j=0; j<nt; ++j) {
            double width = wn*speed[j];
            double strength = factor*boltzmann[(size_t)lower*nt+j]*(-expm1(-c2*wn/t[j]));
            int first = center;
            while (first > 0 && edges[first] > wn-5*width) --first;
            double left = erf((edges[first]-wn)/width);
            for (int i=first; i<bins && edges[i] <= wn+5*width; ++i) {
                double right = erf((edges[i+1]-wn)/width);
                spectrum[(size_t)j*bins+i] += strength*0.5*(right-left);
                left = right;
            }
        }
    }
    fprintf(stderr, "Transitions read: %lu, contributing: %lu\n", count, used);
    if (ferror(stdin) || !count) return 2;
    for (int j=0; j<nt; ++j)
        for (int i=0; i<bins; ++i)
            spectrum[(size_t)j*bins+i] /= edges[i+1]-edges[i];
    return fwrite(spectrum, sizeof(double), (size_t)nt*bins, stdout) != (size_t)nt*bins;
}
