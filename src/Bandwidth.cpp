#include <per4m/Bandwidth.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <fftw3.h>

namespace per4m {
    namespace {

        // ---------- utilities (stats) ----------

        double mean(const std::vector<double> &x) {
            if (x.empty())
                return 0.0;
            long double s = 0.0L;
            for (double v : x)
                s += (long double)v;
            return (double)(s / (long double)x.size());
        }

        double sample_stddev(const std::vector<double> &x) {
            // ddof=1
            const std::size_t n = x.size();
            if (n <= 1)
                return 0.0;
            const double m = mean(x);
            long double acc = 0.0L;
            for (double v : x) {
                long double d = (long double)v - (long double)m;
                acc += d * d;
            }
            acc /= (long double)(n - 1);
            return (double)std::sqrt((double)acc);
        }

        // percentile "linear" (stile type=7), su copia ordinata
        double percentile_linear(std::vector<double> x, double q /*0..100*/) {
            if (x.empty())
                return std::numeric_limits<double>::quiet_NaN();
            if (q <= 0) {
                return *std::min_element(x.begin(), x.end());
            }
            if (q >= 100) {
                return *std::max_element(x.begin(), x.end());
            }
            std::sort(x.begin(), x.end());
            const double p = q / 100.0;
            const double h = (x.size() - 1) * p;
            const auto i = static_cast<std::size_t>(std::floor(h));
            const double frac = h - (double)i;
            if (i + 1 >= x.size())
                return x.back();
            return x[i] + frac * (x[i + 1] - x[i]);
        }

        std::size_t count_unique_exact(std::vector<double> x) {
            if (x.empty())
                return 0;
            std::sort(x.begin(), x.end());
            std::size_t cnt = 1;
            for (std::size_t i = 1; i < x.size(); ++i) {
                if (x[i] != x[i - 1])
                    ++cnt;
            }
            return cnt;
        }

        // ---------- KDEpy-like autogrid ----------
        // Implementazione "molto fedele" come idea:
        // boundary = max(boundary_abs, boundary_rel * range).
        std::vector<double> autogrid_1d(const std::vector<double> &data, double boundary_abs, double boundary_rel, std::size_t num_points) {
            if (data.empty())
                throw std::invalid_argument("autogrid: empty data");
            auto [mn_it, mx_it] = std::minmax_element(data.begin(), data.end());
            const double mn = *mn_it;
            const double mx = *mx_it;
            const double R = mx - mn;
            const double b = std::max(boundary_abs, boundary_rel * R);

            const double left = mn - b;
            const double right = mx + b;

            std::vector<double> grid(num_points);
            if (num_points == 1) {
                grid[0] = 0.5 * (left + right);
                return grid;
            }
            const double dx = (right - left) / (double)(num_points - 1);
            for (std::size_t i = 0; i < num_points; ++i) {
                grid[i] = left + dx * (double)i;
            }
            return grid;
        }

        // ---------- KDEpy-like linear_binning ----------
        // Binning lineare su griglia equispaziata.
        std::vector<double> linear_binning_1d(const std::vector<double> &data, const std::vector<double> &xmesh, const std::vector<double> *weights) {
            const std::size_t n = xmesh.size();
            if (n < 2)
                throw std::invalid_argument("linear_binning: xmesh too small");
            std::vector<double> bins(n, 0.0);

            const double x0 = xmesh.front();
            const double dx = (xmesh.back() - xmesh.front()) / (double)(n - 1);
            if (!(dx > 0.0))
                throw std::invalid_argument("linear_binning: non-positive dx");

            long double wsum = 0.0L;

            for (std::size_t j = 0; j < data.size(); ++j) {
                const double x = data[j];
                const double w = (weights ? (*weights)[j] : 1.0);
                if (w <= 0.0)
                    continue;

                const double u = (x - x0) / dx;
                const long long i = (long long)std::floor(u);

                // outside: KDEpy linear_binning in pratica ignora fuori griglia.
                if (i < 0 || i >= (long long)(n - 1)) {
                    continue;
                }

                const double frac = u - (double)i; // in [0,1)
                const double wl = (1.0 - frac) * w;
                const double wr = frac * w;

                bins[(std::size_t)i] += wl;
                bins[(std::size_t)i + 1] += wr;
                wsum += (long double)w;
            }

            if (wsum <= 0.0L)
                throw std::invalid_argument("linear_binning: no positive weight mass");

            // Normalizza a somma 1, come in KDEpy (assert initial_data.sum()==1)
            for (double &b : bins)
                b = (double)((long double)b / wsum);

            return bins;
        }

        // ---------- ISJ core: fixed point ----------

        long double fixed_point(long double t, long double N, const std::vector<long double> &I_sq, const std::vector<long double> &a2) {
            // fedele al codice Python KDEpy
            const int ell = 7;

            const long double pi = acosl(-1.0L);

            auto sum_term = [&](int power_s, long double time) -> long double {
                // (0.5) * pi^(2s) * sum( I_sq^s * a2 * exp(-I_sq*pi^2*time) )
                long double acc = 0.0L;
                const long double pi2 = pi * pi;
                const long double coeff = 0.5L * powl(pi, (long double)(2 * power_s));
                for (std::size_t i = 0; i < a2.size(); ++i) {
                    // I_sq already squared indices (1..n-1)^2, and stored as long double
                    const long double exp_arg = -I_sq[i] * pi2 * time;
                    const long double e = expl(exp_arg);
                    // I_sq^s
                    const long double isq_pow = powl(I_sq[i], (long double)power_s);
                    acc += isq_pow * a2[i] * e;
                }
                return coeff * acc;
            };

            long double f = sum_term(ell, t);

            if (!(f > 0.0L))
                return -1.0L;

            for (int s = ell - 1; s >= 2; --s) {
                // odd_numbers_prod = prod(1,3,5,...,2s-1,2s+1?) Python: np.arange(1,2*s+1,2)
                long double odd_prod = 1.0L;
                for (int k = 1; k < 2 * s + 1; k += 2)
                    odd_prod *= (long double)k;

                const long double K0 = odd_prod / sqrtl(2.0L * pi);
                const long double cst = (1.0L + powl(0.5L, (long double)(s + 0.5L))) / 3.0L;
                const long double time = powl((2.0L * cst * K0 / (N * f)), (2.0L / (3.0L + 2.0L * (long double)s)));

                f = sum_term(s, time);
                if (!(f > 0.0L))
                    return -1.0L;
            }

            // t_opt = (2*N*sqrt(pi)*f)^(-2/5)
            const long double t_opt = powl(2.0L * N * sqrtl(pi) * f, -2.0L / 5.0L);
            return t - t_opt;
        }

        // ---------- Brent's method (brentq) ----------
        long double brentq(std::function<long double(long double)> f, long double a, long double b, long double xtol = 1e-12L, long double rtol = 1e-10L,
                           int maxiter = 100) {
            long double fa = f(a);
            long double fb = f(b);

            if ((fa > 0.0L && fb > 0.0L) || (fa < 0.0L && fb < 0.0L)) {
                throw std::runtime_error("brentq: f(a) and f(b) must have opposite signs");
            }

            if (fabsl(fa) < fabsl(fb)) {
                std::swap(a, b);
                std::swap(fa, fb);
            }

            long double c = a;
            long double fc = fa;
            bool mflag = true;
            long double d = 0.0L;

            for (int i = 0; i < maxiter; ++i) {
                if (fabsl(b - a) < xtol + rtol * fabsl(b)) {
                    return b;
                }

                if (!(fc == fc) || fabsl(fc) < 1e-14L) {
                    return b;
                }

                long double s;
                if (fabsl(fa - fc) > 1e-14L && fabsl(fb - fc) > 1e-14L) {
                    // inverse quadratic interpolation
                    s = a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc)) + c * fa * fb / ((fc - fa) * (fc - fb));
                } else {
                    // secant method
                    s = b - fb * (b - a) / (fb - fa);
                }

                long double tol_s = xtol + rtol * fabsl(b);
                if ((s <= (3.0L * a + b) / 4.0L || s >= b) || (mflag && fabsl(s - b) >= fabsl(b - c) / 2.0L) ||
                    (!mflag && fabsl(s - b) >= fabsl(c - d) / 2.0L)) {
                    s = (a + b) / 2.0L;
                    mflag = true;
                } else {
                    mflag = false;
                }

                long double fs = f(s);
                d = c;
                c = b;
                fc = fb;

                if (fa * fs < 0.0L) {
                    b = s;
                    fb = fs;
                } else {
                    a = s;
                    fa = fs;
                }

                if (fabsl(fa) < fabsl(fb)) {
                    std::swap(a, b);
                    std::swap(fa, fb);
                }
            }

            throw std::runtime_error("brentq: maximum iterations exceeded");
        }

        // ---------- root finding (faithful strategy) ----------
        // KDEpy usa brentq su [0, tol] e aumenta tol finché trova root.
        // Qui facciamo: cerchiamo un intervallo [a,b] con cambio segno;
        // poi bisection robusta (la brent vera si può aggiungere, ma così è già solidissimo).

        long double root_find_isj(std::size_t N_unique, const std::vector<long double> &I_sq, const std::vector<long double> &a2) {
            // clamp N come in KDEpy
            long double N = (long double)N_unique;
            if (N < 50.0L)
                N = 50.0L;
            if (N > 1050.0L)
                N = 1050.0L;

            long double tol = 1.0e-11L + 0.01L * (N - 50.0L) / 1000.0L;

            auto f = [&](long double t) {
                auto q = fixed_point(t, N_unique, I_sq, a2);
                return q;
            };

            long double a = 0.0L;
            long double fa = f(a);

            // Trova intervallo [a,b] con cambio segno
            long double b = tol;
            long double fb = f(b);

            int attempts = 0;
            while ((fa <= 0.0L && fb <= 0.0L) || (fa >= 0.0L && fb >= 0.0L) || !(fb == fb)) {
                tol *= 2.0L;
                if (tol >= 1.0L) {
                    throw std::runtime_error("ISJ root finding did not converge. Need more data.");
                }
                b = tol;
                fb = f(b);
                if (++attempts > 200) {
                    throw std::runtime_error("ISJ root bracketing failed (too many attempts).");
                }
            }

            // Usa brentq su [a,b]
            long double root = brentq(f, a, b, /*xtol=*/1e-12L, /*rtol=*/1e-10L, /*maxiter=*/100);
            std::cerr << "brentq found root: " << root << " after " << attempts << " attempts to bracket\n";

            if (!(root > 0.0L))
                throw std::runtime_error("ISJ root finding failed to find positive solution.");
            return root;
        }

        // ---------- IQR helpers ----------
        double robust_sigma_min_std_iqr(const std::vector<double> &x) {
            const double sd = sample_stddev(x);
            const double p75 = percentile_linear(x, 75.0);
            const double p25 = percentile_linear(x, 25.0);
            // same constant as KDEpy comment: norm.ppf(.75)-norm.ppf(.25)=1.3489795003921634
            const double iqr_sigma = (p75 - p25) / 1.3489795003921634;
            return std::min(sd, iqr_sigma);
        }

    } // namespace

    // ---------- public API ----------

    double Bandwidth::scott_1d(const std::vector<double> &x) {
        if (x.empty())
            throw std::invalid_argument("scott_1d: empty data");
        const double sigma = robust_sigma_min_std_iqr(x);
        const double n = (double)x.size();
        // dims=1 => n^(-1/(1+4)) = n^(-1/5)
        return sigma * std::pow(n, -1.0 / 5.0);
    }

    double Bandwidth::silverman_1d(const std::vector<double> &x) {
        if (x.empty())
            throw std::invalid_argument("silverman_1d: empty data");
        const std::size_t n = x.size();
        if (n == 1)
            return 1.0;

        double sigma = robust_sigma_min_std_iqr(x);

        if (sigma > 0.0) {
            return sigma * std::pow((double)n * 3.0 / 4.0, -1.0 / 5.0);
        }

        // fallback come KDEpy: usa percentile 99/1
        const double p99 = percentile_linear(x, 99.0);
        const double p01 = percentile_linear(x, 1.0);
        const double iqr2 = (p99 - p01) / 4.6526957480816815; // comment KDEpy
        if (iqr2 > 0.0) {
            return iqr2 * std::pow((double)n * 3.0 / 4.0, -1.0 / 5.0);
        }
        return 1.0;
    }

#ifdef PER4M_USE_FFTW
    // ---------- FFTW DCT-II wrapper ----------
    // SciPy fftpack.dct(x) default: type=2, norm=None (non normalizzata).
    // FFTW: REDFT10 = DCT-II (non normalizzata) con convenzione compatibile.
    std::vector<double> dct_type2_fftw(const std::vector<double> &in) {
        const int n = (int)in.size();
        std::vector<double> out(n);

        // FFTW richiede buffer mutabili
        std::vector<double> tmp_in = in;

        fftw_plan plan = fftw_plan_r2r_1d(n, tmp_in.data(), out.data(),
                                          FFTW_REDFT10, // DCT-II
                                          FFTW_ESTIMATE);

        if (!plan)
            throw std::runtime_error("FFTW plan creation failed for DCT-II.");

        fftw_execute(plan);
        fftw_destroy_plan(plan);

        return out;
    }

    double Bandwidth::isj_1d(const std::vector<double> &x, const std::optional<std::vector<double>>& weights) {
        if (x.empty())
            throw std::invalid_argument("isj_1d: empty data");
        if (weights && weights->size() != x.size()) {
            throw std::invalid_argument("isj_1d: weights size mismatch");
        }

        // KDEpy: n = 2**10
        const std::size_t n_grid = 1u << 10;

        // Filtra weights<=0 (come KDEpy)
        std::vector<double> data;
        std::vector<double> wpos;
        data.reserve(x.size());
        if (weights)
            wpos.reserve(x.size());

        if (weights) {
            for (std::size_t i = 0; i < x.size(); ++i) {
                if ((*weights)[i] > 0.0) {
                    data.push_back(x[i]);
                    wpos.push_back((*weights)[i]);
                }
            }
        } else {
            data = x;
        }

        if (data.empty())
            throw std::invalid_argument("isj_1d: no data after weight filtering");

        // autogrid(data, boundary_abs=6, num_points=n, boundary_rel=0.5)
        const std::vector<double> xmesh = autogrid_1d(data, /*boundary_abs=*/6.0, /*boundary_rel=*/0.5, n_grid);

        // R = max(data) - min(data)
        auto [mn_it, mx_it] = std::minmax_element(data.begin(), data.end());
        const double mn = *mn_it;
        const double mx = *mx_it;
        const double R = mx - mn;
        if (!(R > 0.0)) {
            // dati costanti: ISJ non ha senso; fallback ragionevole
            return 1.0;
        }

        // N = len(unique(data))  (KDEpy)
        const std::size_t N_unique = count_unique_exact(data);
        if (N_unique == 0)
            throw std::runtime_error("isj_1d: N_unique=0");

        // linear_binning su griglia, normalizzato a somma 1
        const std::vector<double> initial_data = linear_binning_1d(data, xmesh, weights ? &wpos : nullptr);

        // DCT-II
        const std::vector<double> a = dct_type2_fftw(initial_data);

        // I_sq = (1..n-1)^2, FLOAT -> long double
        std::vector<long double> I_sq(n_grid - 1);
        for (std::size_t i = 0; i < n_grid - 1; ++i) {
            long double v = (long double)(i + 1);
            I_sq[i] = v * v;
        }

        // a2 = a[1:]^2
        std::vector<long double> a2(n_grid - 1);
        for (std::size_t i = 0; i < n_grid - 1; ++i) {
            long double v = (long double)a[i + 1];
            a2[i] = v * v;
        }

        // solve t_star
        const long double t_star = root_find_isj(N_unique, I_sq, a2);

        // bandwidth = sqrt(t_star) * R
        const long double bw = sqrtl(t_star) * (long double)R;
        return (double)bw;
    }
#endif

} // namespace per4m