// --------------------------------------------------------------------
// OFFSHORE-TO-NEARSHORE WAVE TRANSFORMATION
// (Enhanced Approach Incorporating Wave Obliquity and Refraction)
//
// Overview:
//   This program processes wave data from an input CSV file and computes
//   nearshore wave parameters at a specified depth. It generates two output files:
//     - "output.csv" – Contains the computed nearshore wave parameters.
//     - "report.txt" – Provides descriptive statistics for both the input and
//                      computed variables.
//   The report includes:
//     * The exact command line used to invoke the program.
//     * Detailed descriptive statistics for each variable (count, mean, standard
//       deviation, minimum, maximum, median, and percentiles at 1%, 10%, 25%, 50%,
//       75%, 90%, and 99%).
//     * A table of annual maxima for swh_offshore and swh_local, with a final row
//       indicating the overall maximum values.
//
//   Important Note on Statistics:
//     For the variables alpha_local, swh_local, mwd_local, Ks, Kr, and Hb,
//     the descriptive statistics (count, mean, stddev, min, max, percentiles)
//     are calculated EXCLUDING any time steps where the computed swh_local
//     is zero. This effectively removes waves originating from the land side
//     or those with zero initial offshore height from these specific statistical
//     summaries. Statistics for all other variables include all time steps.
//
//   For directional wave data (mwd_offshore and mwd_local), a hybrid approach is used:
//     - The circular mean and circular standard deviation are computed using the
//       unit-vector method.
//     - The minimum, maximum, median, and quantiles are calculated using ordinary
//       linear statistics on the wrapped angles (in [0,360)).
//     - For mwd_local, the statistics also exclude time steps where swh_local is zero,
//       consistent with the note above.
//
// USAGE:
//   ./transpose input_csv coast_dir depth_d
//
//   Where:
//     input_csv : CSV input file containing at least the following columns:
//                 datetime, swh, mwd, pp1d (additional columns are ignored)
//     coast_dir : Coastline orientation in degrees (clockwise from North)
//     depth_d   : Local depth in meters
//
// EXPECTED CSV INPUT FORMAT (comma-separated):
//     datetime, swh, mwd, pp1d, [additional columns ignored]
//
// OUTPUT CSV FORMAT (comma-separated):
//     datetime,swh_offshore,mwd_offshore,pp1d,L0,L,kh,alpha_offshore,
//     alpha_local,swh_local,mwd_local,Ks,Kr,Hb
//
// Explanation of computed parameters:
//     L0             : Deep-water wavelength, calculated as (g * T²) / (2π)
//     L              : Local wavelength, solved via Newton-Raphson from
//                      L = L0 * tanh((2π * depth_d) / L)
//     kh             : Product of the wave number (k = 2π / L) and local depth (h)
//     alpha_offshore : Signed offshore wave obliquity (crest-to-coast difference),
//                      computed by considering whether the wave approaches from
//                      the "front" (mwd - 90°) or "back" (mwd + 90°) relative to the
//                      coastline.
//     alpha_local    : Local wave angle after refraction, computed from the offshore
//                      obliquity and the refractive factor (tanhl(kh)).
//     mwd_local      : Local mean wave direction, obtained by adjusting the offshore
//                      mwd with the difference between alpha_offshore and alpha_local.
//     Ks             : Shoaling coefficient
//     Kr             : Refraction coefficient
//     Hb             : Breaking wave height (per Miche, 1944), computed as
//                      Hb = 0.142 * L * tanh((2π * depth_d) / L)
//     swh_local      : Local significant wave height, computed as the minimum of
//                      (swh * Ks * Kr) and Hb
//
//   Note: Waves arriving from directions between coast_dir and coast_dir+180°
//         (i.e., from the land side) result in swh_local being set to zero.
//         Initial offshore waves with swh <= 0 also result in swh_local = 0.
//
// Report File Details:
//   The report.txt file includes:
//     - The exact command line used to run the program.
//     - Detailed descriptive statistics for each variable, including count,
//       mean, standard deviation, minimum, maximum, median, and percentiles at
//       1%, 10%, 25%, 50%, 75%, 90%, and 99%. (See note above regarding filtering
//       for specific local wave parameters).
//     - A table displaying the annual maxima for swh_offshore and swh_local,
//       with the final row indicating the overall maximum values.
//
// Compilation Details:
//   To compile the program, use the following command:
//
//     g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic
//         -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++
//         -o transpose transpose.cpp -lm
//
//   Explanation of compile options:
//     - -O3                   : Enables high-level optimizations for maximum performance.
//     - -fopenmp              : Enables OpenMP support for multi-threading.
//     - -march=native         : Optimizes the code for the architecture of the compiling machine.
//     - -std=c++17            : Uses the C++17 standard.
//     - -Wall -Wextra -pedantic: Activates a broad set of compiler warnings to ensure code quality.
//     - -Wconversion          : Warns about implicit type conversions.
//     - -Wsign-conversion     : Warns about implicit sign conversions.
//     - -static, -static-libgcc, -static-libstdc++: Links libraries statically, enhancing portability.
//     - -lm                   : Explicitly link the math library (sometimes needed).
// --------------------------------------------------------------------

#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <cstdlib>
#include <map>
#include <mutex>
#include <set> // For checking variable names efficiently
#include <cstddef> // For size_t

using namespace std;

// --------------------------------------------------------------------
// Physical constants and iteration parameters
// --------------------------------------------------------------------
static const long double G = 9.80665L;               // Acceleration due to gravity (m/s²)
static const long double PI = 3.141592653589793238L; // Mathematical constant π
static const int MAX_ITER = 20;                      // Maximum iterations for Newton-Raphson solver
static const long double TOLERANCE = 1e-12L;         // Convergence tolerance

// --------------------------------------------------------------------
// Utility functions for angle conversion
// --------------------------------------------------------------------
inline long double deg2rad(long double deg)
{
    return deg * (PI / 180.0L);
}

inline long double rad2deg(long double rad)
{
    return rad * (180.0L / PI);
}

// --------------------------------------------------------------------
// Calculate signed offshore obliquity (alpha_offshore in degrees).
//
// For a given offshore mean wave direction (mwd_deg) and a coastline
// direction (coast_deg), we choose whether to compute the crest as
// mwd-90 or mwd+90 depending on the relative orientation. The returned
// value is the signed difference (in degrees) between the wave crest
// and the coastline direction. This sign is later used to adjust the
// rotation due to refraction.
// --------------------------------------------------------------------
long double calcAlphaOffshoreSigned(long double mwd_deg, long double coast_deg)
{
    // Compute relative direction of the offshore wave (0-360)
    long double relativeDir = fmodl((mwd_deg - coast_deg) + 360.0L, 360.0L);
    // Decide which way to compute the crest: use mwd-90 if the wave is coming
    // from the "front" (relativeDir < 180) or mwd+90 if from behind.
    long double crest;
    if (relativeDir < 180.0L) // Wave from land side
        crest = mwd_deg - 90.0L; // Still calculate, but swh_local will be 0 later
    else // Wave from sea side
        crest = mwd_deg + 90.0L;
    // Normalize crest to [0,360)
    crest = fmodl(crest + 360.0L, 360.0L);
    // The signed obliquity is the difference between the crest and the coast.
    long double diff = crest - coast_deg;
    // Normalize difference to [-180, 180)
    while (diff >= 180.0L)
        diff -= 360.0L;
    while (diff < -180.0L)
        diff += 360.0L;
    return diff;
}

// --------------------------------------------------------------------
// Deep-water wavelength calculation.
// --------------------------------------------------------------------
long double deepWaterLength(long double T)
{
    return (T > 0.0L) ? (G * T * T) / (2.0L * PI) : 0.0L;
}

// --------------------------------------------------------------------
// Local wavelength computation using Newton-Raphson.
// Uses an approximation for kh based on k0h for the initial guess.
// --------------------------------------------------------------------
long double localWavelength(long double T, long double depth)
{
    if (T <= 0.0L || depth <= 0.0L)
        return 0.0L;
    const long double L0 = deepWaterLength(T);
    if (L0 <= 0.0L) return 0.0L;

    long double k0h = (2.0L * PI * depth) / L0;
    long double L; // Initial guess for L

    // --- New initial guess based on user formula for kh_approx ---
    // kh_approx ≈ tanh((6/5)^k0h * sqrt(k0h))
    long double kh_approx = k0h / tanhl(powl(6.0L / 5.0L, k0h) * sqrtl(k0h));

    if (kh_approx > TOLERANCE) {
        // If kh = kh_approx, then L = 2*pi / k = 2*pi*depth / kh_approx
        L = (2.0L * PI * depth) / kh_approx;
    } else {
        // Fallback if kh_approx is too small or zero
        // Use Eckart (1952) approximation: L ≈ L0 * sqrt(tanh((k0h)^1.25))
        long double L_eckart = L0 * sqrtl(tanhl(powl(k0h, 1.25L)));
        if (L_eckart > TOLERANCE) {
            L = L_eckart;
        } else {
            L = L0; // Deep water as last resort
        }
    }

    // Wavelength L must be less than or equal to deep water wavelength L0.
    // Cap the initial guess if it exceeds L0.
    if (L > L0) {
        L = L0;
    }
    // Ensure the initial guess L is positive for the iteration.
    if (L <= TOLERANCE) {
        // Use a simpler guess if the primary one failed badly
        L = L0 * tanhl(sqrtl(k0h));
        if (L <= TOLERANCE) L = L0; // Absolute fallback to L0
    }
    // --- End of initial guess calculation ---


    // Newton-Raphson iteration starts here using the initial guess L
    for (int i = 0; i < MAX_ITER; ++i)
    {
        long double current_kh = (2.0L * PI * depth) / L; // Current estimate of kh
        long double tanh_kh = tanhl(current_kh);
        long double F = L - L0 * tanh_kh; // Function F(L) = L - L0*tanh(2*pi*d/L) = 0
        long double sech2_kh = 1.0L - tanh_kh * tanh_kh;
        // Derivative dF/dL = 1 - L0 * d(tanh(kh))/dL
        // d(tanh(kh))/dL = sech^2(kh) * d(kh)/dL = sech^2(kh) * d(2*pi*d/L)/dL
        //                = sech^2(kh) * (-2*pi*d / L^2)
        // dF/dL = 1 - L0 * sech^2(kh) * (-2*pi*d / L^2)
        //       = 1 + L0 * (2*pi*d / L^2) * sech^2(kh)
        long double dF = 1.0L + L0 * (2.0L * PI * depth) / (L * L) * sech2_kh;

        if (fabsl(dF) < TOLERANCE) // Avoid division by zero or near-zero
            break; // Cannot improve further

        long double L_new = L - F / dF; // Newton-Raphson step

        if (L_new <= 0.0L) { // Prevent non-physical results (negative L)
             // If L_new becomes non-positive, the current L might be the best guess,
             // or the iteration might be diverging. Stop iterating.
             break;
        }

        if (fabsl(L_new - L) < TOLERANCE) // Check for convergence
        {
            L = L_new; // Converged
            break;
        }
        L = L_new; // Update L for the next iteration
    }
    // Final check to ensure non-negative return
    return (L < 0.0L ? 0.0L : L);
}

// --------------------------------------------------------------------
// Shoaling coefficient computation.
// --------------------------------------------------------------------
long double shoalingCoefficient(long double k, long double depth)
{
    if (k <= 0.0L || depth <= 0.0L)
        return 1.0L; // Return 1 for invalid inputs or deep water limit
    long double kh = k * depth;
    // Use the standard formula derived from energy flux conservation:
    // Ks = sqrt(Cg0 / Cg) = sqrt( (0.5 * C0) / (n * C) )
    // n = 0.5 * (1 + 2*kh / sinh(2*kh))
    // C/C0 = tanh(kh)
    // Ks = sqrt( 1 / (tanh(kh) * (1 + 2*kh / sinh(2*kh))) )

    long double sinh_2kh = sinhl(2.0L * kh);
    if (fabsl(sinh_2kh) < TOLERANCE) { // Avoid division by zero for very small kh
         // For kh -> 0, Ks^2 ~ 1 / (2*kh). This diverges.
         // Physically, as depth -> 0, waves break. Ks isn't the only factor.
         // Return 1.0 as a safe default if calculation fails near zero.
         return 1.0L;
    }
    long double tanh_kh = tanhl(kh);
    if (fabsl(tanh_kh) < TOLERANCE) return 1.0L; // Avoid division by zero

    long double n = 0.5L * (1.0L + 2.0L * kh / sinh_2kh);
    long double Ks_sq_denominator = tanh_kh * 2.0L * n;

    if (Ks_sq_denominator <= TOLERANCE) return 1.0L; // Avoid division by zero or sqrt of negative

    long double Ks_sq = 1.0L / Ks_sq_denominator;

    return (Ks_sq > 0.0L) ? sqrtl(Ks_sq) : 1.0L; // Return 1 if calculation yields non-positive result
}

// --------------------------------------------------------------------
// Structure for descriptive statistics
// --------------------------------------------------------------------
struct DescriptiveStats
{
    size_t count;
    long double mean;    // For directional data, circular mean will be used
    long double stddev;  // For directional data, circular std. dev. will be used
    long double min;
    long double p1;     // 1st percentile
    long double p10;    // 10th percentile
    long double p25;    // 25th percentile
    long double median; // 50th percentile
    long double p75;    // 75th percentile
    long double p90;    // 90th percentile
    long double p99;    // 99th percentile
    long double max;
};

// --------------------------------------------------------------------
// Compute descriptive statistics (for linear data).
// --------------------------------------------------------------------
DescriptiveStats computeStats(const vector<long double> &data)
{
    DescriptiveStats stats = {0, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L};
    if (data.empty()) {
        stats.min = numeric_limits<long double>::quiet_NaN();
        stats.max = numeric_limits<long double>::quiet_NaN();
        stats.mean = numeric_limits<long double>::quiet_NaN();
        stats.stddev = numeric_limits<long double>::quiet_NaN();
        stats.p1 = numeric_limits<long double>::quiet_NaN();
        stats.p10 = numeric_limits<long double>::quiet_NaN();
        stats.p25 = numeric_limits<long double>::quiet_NaN();
        stats.median = numeric_limits<long double>::quiet_NaN();
        stats.p75 = numeric_limits<long double>::quiet_NaN();
        stats.p90 = numeric_limits<long double>::quiet_NaN();
        stats.p99 = numeric_limits<long double>::quiet_NaN();
        return stats;
    }
    stats.count = data.size();
    long double sum = 0.0L;
    stats.min = data[0];
    stats.max = data[0];
    for (long double x : data)
    {
        sum += x;
        if (x < stats.min)
            stats.min = x;
        if (x > stats.max)
            stats.max = x;
    }
    stats.mean = sum / stats.count;
    long double variance = 0.0L;
    for (long double x : data)
    {
        long double diff = x - stats.mean;
        variance += diff * diff;
    }
    stats.stddev = (stats.count > 1) ? sqrtl(variance / (stats.count - 1)) : 0.0L;
    vector<long double> sortedData = data;
    sort(sortedData.begin(), sortedData.end());
    auto getPercentile = [&](long double p) -> long double {
        if (sortedData.empty()) return numeric_limits<long double>::quiet_NaN();
        long double pos = p * (sortedData.size() - 1);
        size_t idx = static_cast<size_t>(floorl(pos));
        long double frac = pos - idx;
        if (idx + 1 < sortedData.size()) {
             return sortedData[idx] * (1.0L - frac) + sortedData[idx + 1] * frac;
        } else {
             return sortedData[idx]; // Handles case where size=1 or p=1.0
        }
    };
    stats.p1     = getPercentile(0.01L);
    stats.p10    = getPercentile(0.10L);
    stats.p25    = getPercentile(0.25L);
    stats.median = getPercentile(0.50L);
    stats.p75    = getPercentile(0.75L);
    stats.p90    = getPercentile(0.90L);
    stats.p99    = getPercentile(0.99L);
    // Ensure min/max are consistent if only one element
    if (stats.count == 1) {
        stats.min = sortedData[0];
        stats.max = sortedData[0];
    }
    return stats;
}

// --------------------------------------------------------------------
// Compute hybrid circular statistics for directional data.
// --------------------------------------------------------------------
DescriptiveStats computeHybridCircularStats(const vector<long double> &data)
{
    if (data.empty()) {
        DescriptiveStats stats = {0, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L};
        stats.min = numeric_limits<long double>::quiet_NaN();
        stats.max = numeric_limits<long double>::quiet_NaN();
        stats.mean = numeric_limits<long double>::quiet_NaN(); // Circular mean is NaN
        stats.stddev = numeric_limits<long double>::quiet_NaN();// Circular stddev is NaN
        stats.p1 = numeric_limits<long double>::quiet_NaN();
        stats.p10 = numeric_limits<long double>::quiet_NaN();
        stats.p25 = numeric_limits<long double>::quiet_NaN();
        stats.median = numeric_limits<long double>::quiet_NaN();
        stats.p75 = numeric_limits<long double>::quiet_NaN();
        stats.p90 = numeric_limits<long double>::quiet_NaN();
        stats.p99 = numeric_limits<long double>::quiet_NaN();
        return stats;
    }

    vector<long double> wrapped;
    wrapped.reserve(data.size());
    for (long double d : data) {
        long double w = fmodl(d, 360.0L);
        if (w < 0)
            w += 360.0L;
        wrapped.push_back(w);
    }
    // Calculate linear stats (min, max, percentiles) on wrapped data first
    DescriptiveStats hybrid = computeStats(wrapped); // This sets count, min, max, percentiles

    // Calculate circular mean and stddev
    long double sumSin = 0.0L, sumCos = 0.0L;
    for (long double d_wrapped : wrapped) { // Use the wrapped data for consistency
        long double rad = deg2rad(d_wrapped);
        sumSin += sinl(rad);
        sumCos += cosl(rad);
    }

    long double meanRad = atan2l(sumSin, sumCos); // Result is in [-PI, PI]
    if (meanRad < 0)
        meanRad += 2.0L * PI; // Shift to [0, 2*PI)
    long double circMean = rad2deg(meanRad); // Convert to degrees [0, 360)

    long double R_bar = sqrtl(sumSin * sumSin + sumCos * sumCos) / wrapped.size(); // Mean resultant length
    long double circStdRad = sqrtl(-2.0L * logl(max(R_bar, static_cast<long double>(1e-15L)))); // Use max to avoid log(0) if R_bar is exactly 0
    long double circStd = rad2deg(circStdRad); // Circular standard deviation in degrees

    // Update the stats structure with circular values
    hybrid.mean = circMean;
    hybrid.stddev = circStd;

    return hybrid;
}


// --------------------------------------------------------------------
// Main function
// --------------------------------------------------------------------
int main(int argc, char *argv[])
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout << fixed << setprecision(10); // Apply precision globally for cout if needed
    cerr << fixed << setprecision(4);  // Apply precision globally for cerr if needed

    if (argc != 4)
    {
        cerr << "Usage: " << argv[0] << " input_csv coast_dir depth_d\n"
             << "Example: " << argv[0] << " input.csv 10 5\n";
        return 1;
    }
    string input_csv = argv[1];
    long double coast_dir_input, depth_d;
    try {
        coast_dir_input = stold(argv[2]);
        depth_d = stold(argv[3]);
    }
    catch (...) {
        cerr << "Error: cannot parse coast_dir or depth_d as numbers.\n";
        return 1;
    }
     // Normalize coast_dir to [0, 360)
    long double coast_dir = fmodl(coast_dir_input, 360.0L);
    if (coast_dir < 0) coast_dir += 360.0L;

    if (depth_d <= 0.0L) {
        cerr << "Invalid depth value: " << depth_d << ". Must be positive.\n";
        return 1;
    }
    ifstream inFile(input_csv);
    if (!inFile.is_open()) {
        cerr << "ERROR: unable to open input file " << input_csv << "\n";
        return 1;
    }
    ofstream outFile("output.csv");
    if (!outFile.is_open()) {
        cerr << "ERROR: unable to create output.csv\n";
        return 1;
    }
    // Apply precision to output file stream
    outFile << fixed << setprecision(10);

    // Write output CSV header
    outFile << "datetime,"
            << "swh_offshore,mwd_offshore,pp1d,"
            << "L0,L,kh,alpha_offshore,alpha_local,"
            << "swh_local,mwd_local,Ks,Kr,Hb\n";

    // Preallocate arrays for statistics (13 columns)
    const size_t NUM_COLS = 13;
    vector<vector<long double>> statsData(NUM_COLS); // No initial size, will resize later

    // Read and discard CSV header
    string header;
    getline(inFile, header);
    // Read remaining lines into a vector
    vector<string> lines;
    string line;
    while(getline(inFile, line)) {
        // Basic check for empty or whitespace-only lines
        if (!line.empty() && line.find_first_not_of(" \t\n\v\f\r") != string::npos)
            lines.push_back(line);
    }
    inFile.close();

    if (lines.empty()) {
        cerr << "Warning: Input file contains no data lines after the header.\n";
        outFile.close(); // Close the output file (it just has the header)
        // Create an empty report file
        ofstream reportFile("report.txt");
        if (!reportFile.is_open()) {
             cerr << "ERROR: unable to create report.txt\n";
             return 1; // Still an error if report can't be made
        }
              
        // Write the exact command line used to invoke the program
        reportFile << "Command line:";
        for (int i = 0; i < argc; ++i) { reportFile << " " << argv[i]; }
        reportFile << "\n\n"; // Add extra newline for spacing before the next section
    
        reportFile << "Input file contained no valid data.\n";
        reportFile.close();
        return 0; // Not an error, just no data to process
    }


    // Sort lines by datetime (first token) and remove duplicates based on datetime.
    sort(lines.begin(), lines.end(), [](const string &a, const string &b) {
        size_t posA = a.find(',');
        size_t posB = b.find(',');
        string dtA = (posA == string::npos) ? a : a.substr(0, posA);
        string dtB = (posB == string::npos) ? b : b.substr(0, posB);
        return dtA < dtB;
    });
    vector<string> uniqueLines;
    uniqueLines.reserve(lines.size());
    if (!lines.empty()) {
        uniqueLines.push_back(lines[0]);
        for (size_t i = 1; i < lines.size(); ++i) {
            size_t pos = lines[i].find(',');
            string dt = (pos == string::npos) ? lines[i] : lines[i].substr(0, pos);
            size_t posLast = uniqueLines.back().find(',');
            string lastDt = (posLast == string::npos) ? uniqueLines.back() : uniqueLines.back().substr(0, posLast);
            if (lastDt != dt) {
                uniqueLines.push_back(lines[i]);
            }
            // else: skip duplicate datetime
        }
    }
    lines.clear(); // Free memory from original potentially larger vector

    // Resize statsData based on the number of unique lines
    size_t num_records = uniqueLines.size();
    for (size_t j = 0; j < NUM_COLS; j++)
        statsData[j].resize(num_records, 0.0L);

    // Preallocate output storage.
    vector<string> outputLines(num_records);


    // Prepare thread-local annual maximum maps.
    size_t numThreads = static_cast<size_t>(omp_get_max_threads());
    vector< map<string, long double> > localAnnualMaxOffshore(numThreads);
    vector< map<string, long double> > localAnnualMaxLocal(numThreads);

    // ----------------------------------------------------------------
    // Process the unique lines in parallel.
    // ----------------------------------------------------------------
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_records; ++i) {
        size_t tid = static_cast<size_t>(omp_get_thread_num());
        const auto &l = uniqueLines[i];
        stringstream ss(l);
        string segment;
        vector<string> fields;
        while (getline(ss, segment, ',')) {
            fields.push_back(segment);
        }

        string datetime = "";
        long double swh = 0.0L, mwd = 0.0L, pp1d = 0.0L;
        bool parse_ok = false;

        if (fields.size() >= 4) {
            try {
                datetime = fields[0];
                swh = stold(fields[1]);
                mwd = stold(fields[2]);
                pp1d = stold(fields[3]);
                parse_ok = true;
            } catch (const std::invalid_argument& ia) {
                // cerr << "Warning: Skipping line " << i+1 << " due to invalid number format: " << ia.what() << "\n";
                parse_ok = false;
            } catch (const std::out_of_range& oor) {
                // cerr << "Warning: Skipping line " << i+1 << " due to out of range number: " << oor.what() << "\n";
                parse_ok = false;
            } catch (...) {
                // cerr << "Warning: Skipping line " << i+1 << " due to unknown parsing error.\n";
                parse_ok = false;
            }
        } else {
             // cerr << "Warning: Skipping line " << i+1 << " due to insufficient columns.\n";
             parse_ok = false;
        }

        vector<long double> record(NUM_COLS, 0.0L); // Initialize with zeros
        ostringstream oss;
        oss << fixed << setprecision(10); // Ensure precision in the output string stream

        if (!parse_ok || swh <= 0.0L || pp1d <= 0.0L) {
            // Handle invalid input or non-physical waves
            oss << datetime << "," << swh << "," << mwd << "," << pp1d;
            oss << ",0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0"; // Zeroes for calculated fields
            record = {swh, mwd, pp1d, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L};
            // Update annual max even if input is zero/invalid (max will be >= 0)
             if (!datetime.empty() && datetime.size() >= 4) {
                 string year = datetime.substr(0, 4);
                 // Use mutex or critical section if map access wasn't thread-local
                 localAnnualMaxOffshore[tid][year] = max(localAnnualMaxOffshore[tid].count(year) ? localAnnualMaxOffshore[tid][year] : -1.0L , swh);
                 localAnnualMaxLocal[tid][year] = max(localAnnualMaxLocal[tid].count(year) ? localAnnualMaxLocal[tid][year] : -1.0L, 0.0L); // swh_local is 0 here
             }
        } else {
            // Valid offshore wave data, proceed with calculation
            long double mwd_offshore_norm = fmodl(mwd, 360.0L);
             if (mwd_offshore_norm < 0) mwd_offshore_norm += 360.0L;

            long double relativeDir = fmodl((mwd_offshore_norm - coast_dir) + 360.0L, 360.0L);

            // Check if wave comes from land side (relative angle between 0 and 180 degrees exclusive of 180)
            // Note: relativeDir = 180 means wave is exactly parallel to coast from seaward side
            if (relativeDir > 0.0L && relativeDir < 180.0L) {
                // Wave from land side: set local parameters to zero
                oss << datetime << "," << swh << "," << mwd_offshore_norm << "," << pp1d;
                oss << ",0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0"; // Zeroes for calculated fields
                record = {swh, mwd_offshore_norm, pp1d, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L};
                 if (!datetime.empty() && datetime.size() >= 4) {
                     string year = datetime.substr(0, 4);
                     localAnnualMaxOffshore[tid][year] = max(localAnnualMaxOffshore[tid].count(year) ? localAnnualMaxOffshore[tid][year] : -1.0L, swh);
                     localAnnualMaxLocal[tid][year] = max(localAnnualMaxLocal[tid].count(year) ? localAnnualMaxLocal[tid][year] : -1.0L, 0.0L); // swh_local is 0
                 }
            } else {
                // Wave from sea side (relative angle 0, or 180 to 360)
                long double L0 = deepWaterLength(pp1d);
                long double L = localWavelength(pp1d, depth_d);
                long double kLocal = (L > TOLERANCE) ? (2.0L * PI / L) : 0.0L;
                long double kh = kLocal * depth_d;

                long double alpha_offshore = calcAlphaOffshoreSigned(mwd_offshore_norm, coast_dir);
                long double alpha_offshore_rad = deg2rad(alpha_offshore);

                // Refraction calculation using Snell's Law: sin(alpha_local)/C = sin(alpha_offshore)/C0
                // C = L/T, C0 = L0/T => C/C0 = L/L0 = tanh(kh)
                // sin(alpha_local) = sin(alpha_offshore) * (C/C0) = sin(alpha_offshore) * tanh(kh)
                long double tanh_kh_val = (kh > 0.0L) ? tanhl(kh) : 0.0L; // Avoid tanh(0) issues if kh is exactly 0
                long double sinAlphaLocal = sinl(alpha_offshore_rad) * tanh_kh_val;


                // Clamp sinAlphaLocal to [-1, 1] to avoid domain errors with asinl
                if (sinAlphaLocal > 1.0L) sinAlphaLocal = 1.0L;
                else if (sinAlphaLocal < -1.0L) sinAlphaLocal = -1.0L;

                long double alpha_local_rad = asinl(sinAlphaLocal);
                long double alpha_local = rad2deg(alpha_local_rad);

                // Adjust local MWD: mwd_local = mwd_offshore - (change in angle relative to normal)
                // Change in angle relative to normal = alpha_offshore - alpha_local
                long double mwd_local = fmodl(mwd_offshore_norm - (alpha_offshore - alpha_local) + 360.0L, 360.0L);

                long double Ks = shoalingCoefficient(kLocal, depth_d);

                // Refraction coefficient Kr = sqrt(cos(alpha_offshore) / cos(alpha_local))
                long double cosAlphaOffshore = cosl(alpha_offshore_rad);
                long double cosAlphaLocal = cosl(alpha_local_rad);
                long double Kr = 1.0L; // Default value
                 // Ensure cosines are non-negative (should be true for angles between -90 and 90)
                 // Avoid division by zero if alpha_local is +/- 90 degrees
                 if (cosAlphaOffshore >= 0 && cosAlphaLocal > TOLERANCE) {
                     Kr = sqrtl(cosAlphaOffshore / cosAlphaLocal);
                 } else if (cosAlphaLocal <= TOLERANCE && cosAlphaOffshore > TOLERANCE) {
                     // Grazing incidence or angle becoming 90 deg. Kr -> infinity?
                     // Physically, energy spreads out. Some models limit Kr or use alternatives.
                     // Setting Kr=1.0 is a simple fallback, but might underestimate focusing.
                     // Let's keep Kr=1.0 as a safe fallback here.
                     Kr = 1.0L; // Or consider a large finite value?
                 }
                 // If cosAlphaOffshore < 0 (unlikely with signed alpha), Kr is complex. Fallback Kr=1.0.


                // Breaking wave height (Miche, 1944) - Hb = 0.142 * L * tanh(kh)
                long double Hb = (L > TOLERANCE && kh > 0.0L) ? 0.142L * L * tanh_kh_val : 0.0L;

                long double swh_transformed = swh * Ks * Kr;
                long double swh_local = (Hb > 0.0L) ? min(swh_transformed, Hb) : swh_transformed; // Apply breaking limit

                // Ensure swh_local is non-negative
                if (swh_local < 0.0L) swh_local = 0.0L;

                oss << datetime << "," << swh << "," << mwd_offshore_norm << "," << pp1d << ","
                    << L0 << "," << L << "," << kh << ","
                    << alpha_offshore << "," << alpha_local << ","
                    << swh_local << "," << mwd_local << ","
                    << Ks << "," << Kr << "," << Hb;

                record = {swh, mwd_offshore_norm, pp1d, L0, L, kh, alpha_offshore,
                          alpha_local, swh_local, mwd_local, Ks, Kr, Hb};

                if (!datetime.empty() && datetime.size() >= 4) {
                    string year = datetime.substr(0, 4);
                    localAnnualMaxOffshore[tid][year] = max(localAnnualMaxOffshore[tid].count(year) ? localAnnualMaxOffshore[tid][year] : -1.0L, swh);
                    localAnnualMaxLocal[tid][year] = max(localAnnualMaxLocal[tid].count(year) ? localAnnualMaxLocal[tid][year] : -1.0L, swh_local);
                }
            }
        }
        outputLines[i] = oss.str();
        // Store record data into the main statsData arrays
        for (size_t j = 0; j < NUM_COLS; j++) {
             // Ensure we don't write out of bounds if record is somehow smaller (shouldn't happen with init)
             if (j < record.size()) {
                 statsData[j][i] = record[j];
             }
        }
    } // end parallel for

    // Merge thread-local annual maximum maps.
    map<string, long double> annualMaxOffshore, annualMaxLocal;
    for (size_t t = 0; t < numThreads; t++) {
        for (const auto &p : localAnnualMaxOffshore[t]) {
            annualMaxOffshore[p.first] = max(annualMaxOffshore.count(p.first) ? annualMaxOffshore[p.first] : -1.0L, p.second);
        }
        for (const auto &p : localAnnualMaxLocal[t]) {
            annualMaxLocal[p.first] = max(annualMaxLocal.count(p.first) ? annualMaxLocal[p.first] : -1.0L, p.second);
        }
    }

    // Write output lines sequentially to ensure order.
    for (const auto &s : outputLines)
        outFile << s << "\n";
    outFile.close();

    // ----------------------------------------------------------------
    // Create report.txt with descriptive statistics and annual maxima.
    // ----------------------------------------------------------------
    const int LINE_WIDTH = 105; // Adjusted slightly for potentially 3 columns
    const int VARIABLES_PER_ROW = 3;
    int colWidth = LINE_WIDTH / VARIABLES_PER_ROW;
    ofstream reportFile("report.txt");
    if (!reportFile.is_open()) {
        cerr << "ERROR: unable to create report.txt\n";
        return 1;
    }
    // Apply precision to report file stream
    reportFile << fixed << setprecision(10);

    reportFile << "Command line: " << argv[0] << " " << argv[1] << " " << argv[2] << " " << argv[3] << "\n\n";
    reportFile << "Descriptive Statistics Report\n";
    reportFile << "(Note: Stats for alpha_local, swh_local, mwd_local, Ks, Kr, Hb exclude steps where swh_local = 0)\n";
    reportFile << string(LINE_WIDTH, '=') << "\n\n";

    vector<string> varNames = {
        "swh_offshore", "mwd_offshore", "pp1d", "L0", "L", "kh",
        "alpha_offshore", "alpha_local", "swh_local", "mwd_local", "Ks", "Kr", "Hb"
    };
    // Set of variables that require filtering (based on swh_local > 0)
    set<string> vars_to_filter = {"alpha_local", "swh_local", "mwd_local", "Ks", "Kr", "Hb"};

    vector<string> statLabels = {"Count", "Mean", "Std. Dev.", "Min", "P1", "P10", "P25",
                                 "Median", "P75", "P90", "P99", "Max"};

    // Find the index of swh_local (used for filtering)
    size_t swh_local_idx = 8; // Based on varNames order

    for (size_t i = 0; i < NUM_COLS; i += VARIABLES_PER_ROW) {
        // Print variable names header row
        for (size_t j = 0; j < VARIABLES_PER_ROW && (i + j) < NUM_COLS; ++j) {
            reportFile << left << setw(colWidth) << varNames[i + j];
        }
        reportFile << "\n";
        reportFile << string(static_cast<size_t>(LINE_WIDTH), '-') << "\n"; // Separator below names

        vector<DescriptiveStats> groupStats;
        vector< vector<long double> > groupData(VARIABLES_PER_ROW); // To hold data (filtered or full) for the group

        // Prepare data (filter if necessary) and calculate stats for the current group of variables
        for (size_t j = 0; j < VARIABLES_PER_ROW && (i + j) < NUM_COLS; ++j) {
            size_t current_col_index = i + j;
            string current_var_name = varNames[current_col_index];

            bool needs_filtering = (vars_to_filter.count(current_var_name) > 0);

            if (needs_filtering) {
                vector<long double> filtered_data;
                filtered_data.reserve(num_records);
                const vector<long double>& swh_local_data = statsData[swh_local_idx];
                const vector<long double>& current_col_data = statsData[current_col_index];

                for(size_t k = 0; k < num_records; ++k) {
                    // Apply the filter condition: swh_local > 0
                    if (swh_local_data[k] > TOLERANCE) { // Use tolerance for float comparison
                        filtered_data.push_back(current_col_data[k]);
                    }
                }
                groupData[j] = filtered_data; // Store filtered data

                // Compute stats on the filtered data
                if (current_var_name == "mwd_local") {
                     groupStats.push_back(computeHybridCircularStats(filtered_data));
                } else {
                     groupStats.push_back(computeStats(filtered_data));
                }
            } else {
                // Use the original (unfiltered) data
                groupData[j] = statsData[current_col_index]; // Store full data
                if (current_var_name == "mwd_offshore") {
                     groupStats.push_back(computeHybridCircularStats(statsData[current_col_index]));
                } else {
                     groupStats.push_back(computeStats(statsData[current_col_index]));
                }
            }
        }

        // Print the statistics for the current group
        for (const auto &label : statLabels) {
            for (size_t j = 0; j < groupStats.size(); ++j) { // groupStats.size() will be <= VARIABLES_PER_ROW
                ostringstream oss;
                oss << fixed << setprecision(10); // Ensure precision in report numbers
                long double value = 0.0;
                bool is_count = false;

                if (label == "Count") {
                    is_count = true;
                } else if (label == "Mean")       value = groupStats[j].mean;
                else if (label == "Std. Dev.")  value = groupStats[j].stddev;
                else if (label == "Min")        value = groupStats[j].min;
                else if (label == "P1")         value = groupStats[j].p1;
                else if (label == "P10")        value = groupStats[j].p10;
                else if (label == "P25")        value = groupStats[j].p25;
                else if (label == "Median")     value = groupStats[j].median;
                else if (label == "P75")        value = groupStats[j].p75;
                else if (label == "P90")        value = groupStats[j].p90;
                else if (label == "P99")        value = groupStats[j].p99;
                else if (label == "Max")        value = groupStats[j].max;

                oss << label << ": ";
                if (is_count) {
                    oss << groupStats[j].count;
                } else {
                   if (isnan(value)) {
                       oss << "NaN";
                   } else {
                       // Adjust precision for specific stats if needed, e.g., angles
                       if (varNames[i+j] == "mwd_offshore" || varNames[i+j] == "mwd_local" || varNames[i+j] == "alpha_offshore" || varNames[i+j] == "alpha_local") {
                           oss << fixed << setprecision(3) << value;
                       } else {
                           oss << fixed << setprecision(6) << value; // Default precision for others
                       }
                   }
                }

                reportFile << left << setw(colWidth) << oss.str();
            }
            reportFile << "\n";
        }
        reportFile << string(static_cast<size_t>(LINE_WIDTH), '-') << "\n\n"; // Separator after each group
    }

    // Annual Maxima table.
    reportFile << "\nAnnual Maxima of swh_offshore and swh_local:\n";
    int tableColWidth = 30;
    reportFile << left << setw(tableColWidth) << "Year"
               << left << setw(tableColWidth) << "swh_offshore (Max)"
               << left << setw(tableColWidth) << "swh_local (Max)" << "\n";
    // Explicit cast to size_t to avoid sign conversion warning
    reportFile << string(static_cast<size_t>(tableColWidth * 3), '-') << "\n";

    long double overallOffshore = -1.0L, overallLocal = -1.0L; // Initialize to -1 to correctly find max >= 0
    vector<string> years;
    for(const auto& pair : annualMaxOffshore) years.push_back(pair.first);
    sort(years.begin(), years.end()); // Sort years for ordered output

    reportFile << fixed << setprecision(6); // Set precision for maxima table

    for (const auto &year : years) {
        long double offVal = annualMaxOffshore[year];
        long double localVal = (annualMaxLocal.count(year)) ? annualMaxLocal[year] : 0.0L; // Default to 0 if no entry for local

        reportFile << left << setw(tableColWidth) << year
                   << left << setw(tableColWidth) << offVal
                   << left << setw(tableColWidth) << localVal << "\n";

        overallOffshore = max(overallOffshore, offVal);
        overallLocal = max(overallLocal, localVal);
    }
     // Explicit cast to size_t to avoid sign conversion warning
     reportFile << string(static_cast<size_t>(tableColWidth * 3), '-') << "\n";
    reportFile << left << setw(tableColWidth) << "Overall Max"
               << left << setw(tableColWidth) << overallOffshore
               << left << setw(tableColWidth) << overallLocal << "\n";

    reportFile.close();
    cout << "Processing complete. Output written to output.csv and report.txt\n";
    return 0;
}