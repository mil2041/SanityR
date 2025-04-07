#include "PsiGamma.h"
#include <cmath>

// Anonymous namespace
namespace {

// Constants
constexpr double EULER_MASCHERONI = -0.577215664901532;
constexpr double PI_SQUARED_OVER_6 = 1.644934066848226;

// Helper function for series expansion terms
inline double series_sum(double x, const double* coefficients, int n_terms, int start_power) {
    double sum = 0.0;
    double x_i = std::pow(x, start_power);
    double x_square = x * x;

    for (int i = 0; i < n_terms; ++i) {
        sum += coefficients[i] / x_i;
        x_i *= x_square;  // Skip every other power
    }
    return sum;
}

}

namespace psi {

double Psi_0(double x) {
    if (x == 1.0) return EULER_MASCHERONI;

    // Coefficients for digamma series expansion
    static constexpr double coeffs[] = {
        -1.0/12,      // 1/x^2
         1.0/120,     // 1/x^4
        -1.0/252,     // 1/x^6
         1.0/240,     // 1/x^8
        -5.0/660,     // 1/x^10
         691.0/32760, // 1/x^12
        -1.0/12       // 1/x^14
    };

    const double log_x = std::log(x);
    const double series = series_sum(x, coeffs, 7, 2);

    return log_x - 0.5/x + series;
}

double Psi_1(double x) {
    if (x == 1.0) return PI_SQUARED_OVER_6;

    // Coefficients for trigamma series expansion
    static constexpr double coeffs[] = {
         1.0,        // 1/x
         1.0/6,      // 1/x^3
        -1.0/30,     // 1/x^5
         1.0/42,     // 1/x^7
        -1.0/30,     // 1/x^9
        -5.0/66,     // 1/x^11
        -691.0/2730, // 1/x^13
         7.0/6       // 1/x^15
    };

    const double inv_x2 = 0.5 / (x*x);
    const double series = series_sum(x, coeffs, 8, 1);

    return inv_x2 + series;
}

} // namespace psi
