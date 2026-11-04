// Native (no Python) throughput benchmark of the gkpbp decoder on a problem instance exported by scripts/export_instance.py.
// Usage: bench_native instance.bin [threads=1] [repeats=3]
// Binary layout (little endian): int32 m,n,E,B,max_iter,osd_order ; float64 ms_scaling ; int32 indptr[m+1], indices[E] ;
//                                uint8 syndromes[B*m] ; float64 llr[B*n]  (per-shot priors)
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "gkpbp/decoder.hpp"

template <class T> static void rd(FILE* f, T* p, size_t k) { if (fread(p, sizeof(T), k, f) != k) { fprintf(stderr, "short read\n"); exit(2); } }

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s instance.bin [threads] [repeats]\n", argv[0]); return 1; }
    const int threads = argc > 2 ? atoi(argv[2]) : 1, repeats = argc > 3 ? atoi(argv[3]) : 3;
    FILE* f = fopen(argv[1], "rb"); if (!f) { perror("open"); return 1; }
    int32_t hdr[6]; double ms; rd(f, hdr, 6); rd(f, &ms, 1);
    const int m = hdr[0], n = hdr[1], E = hdr[2], B = hdr[3], max_iter = hdr[4], osd = hdr[5];
    std::vector<int32_t> indptr(m + 1), indices(E); std::vector<uint8_t> S((size_t)B * m); std::vector<double> L((size_t)B * n);
    rd(f, indptr.data(), m + 1); rd(f, indices.data(), E); rd(f, S.data(), S.size()); rd(f, L.data(), L.size()); fclose(f);
    std::vector<uint8_t> out((size_t)B * n); std::vector<char> conv(B), ok(B); std::vector<int> iters(B);
    double best = 1e30;
    for (int rep = 0; rep < repeats; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
        {
            gkpbp::BpOsd dec(m, n, indptr.data(), indices.data(), max_iter, ms, osd);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 4)
#endif
            for (int b = 0; b < B; ++b) {
                auto r = dec.decode(&S[(size_t)b * m], &L[(size_t)b * n], &out[(size_t)b * n]);
                conv[b] = r.bp_converged; ok[b] = r.syndrome_ok; iters[b] = r.iterations;
            }
        }
        best = std::min(best, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    }
    double c = 0, o = 0, it = 0; for (int b = 0; b < B; ++b) { c += conv[b]; o += ok[b]; it += iters[b]; }
    printf("m=%d n=%d E=%d B=%d osd_order=%d threads=%d | best of %d: %.4f s  %.1f us/shot  %.0f shots/s | BP converged %.3f  syndrome ok %.3f  mean iters %.1f\n",
           m, n, E, B, osd, threads, repeats, best, best / B * 1e6, B / best, c / B, o / B, it / B);
    return 0;
}
