// dllmain.cpp : Defines the entry point for the DLL application.
#include "pch.h"
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>

const int versionMajor = 0;
const int versionMinor = 0;
const int versionFix = 1;

#define METHOD_EXPORTS
#ifdef METHOD_EXPORTS
#define EXPORT __declspec(dllexport)
#else
#define EXPORT __declspec(dllimport)
#endif

const int MASS_RANGE = 5000;                                // Encoding values up to 5000 m/z
const int MASS_MULTIPLIER = 100;                            // Encoding values with 0.01 precision
const int ENCODING_SIZE = MASS_RANGE * MASS_MULTIPLIER;     // The total length of an encoding vector
const int APPROX_NNZ_PER_ROW = 100;                         // Approximate number of ions assumed
const int ROUNDING_ACCURACY = 1000;                         // Rounding precision for converting f32 to i32, the exact precision is (int) round(val * 1000.0f)
const double ONE_OVER_SQRT_PI = 0.39894228040143267793994605993438;

extern "C" {
    EXPORT int* findTopCandidates(int*, int*, 
                                  int*, int*, 
                                  int, int, 
                                  int, int, 
                                  int, float,
                                  bool, bool,
                                  int, int);

    EXPORT int* findTopCandidatesInt(int*, int*,
                                     int*, int*,
                                     int, int,
                                     int, int,
                                     int, float,
                                     bool, bool,
                                     int, int);

    EXPORT int* findTopCandidates2(int*, int*,
                                   int*, int*,
                                   int, int,
                                   int, int,
                                   int, float,
                                   bool, bool,
                                   int, int);

    EXPORT int* findTopCandidates2Int(int*, int*,
                                      int*, int*,
                                      int, int,
                                      int, int,
                                      int, float,
                                      bool, bool,
                                      int, int);

    EXPORT int* findTopCandidatesBatched(int*, int*,
                                         int*, int*,
                                         int, int,
                                         int, int,
                                         int, float,
                                         bool, bool,
                                         int,
                                         int, int);

    EXPORT int* findTopCandidatesBatchedInt(int*, int*,
                                            int*, int*,
                                            int, int,
                                            int, int,
                                            int, float,
                                            bool, bool,
                                            int,
                                            int, int);

    EXPORT int* findTopCandidatesBatched2(int*, int*,
                                          int*, int*,
                                          int, int,
                                          int, int,
                                          int, float,
                                          bool, bool,
                                          int,
                                          int, int);

    EXPORT int* findTopCandidatesBatched2Int(int*, int*,
                                             int*, int*,
                                             int, int,
                                             int, int,
                                             int, float,
                                             bool, bool,
                                             int,
                                             int, int);

    EXPORT int releaseMemory(int*);
}

float squared(float);
float normpdf(float, float, float);

/// <summary>
/// Exception thrown for functions that haven't been implemented yet.
/// </summary>
class NotImplementedException : public std::logic_error
{
public:
    NotImplementedException() : std::logic_error("Function not yet implemented!") { };
};

/// <summary>
/// A function that calculates the top n candidates for each spectrum. Maps to CPU_SVf32 (and variations).
/// </summary>
/// <param name="candidatesValues">An integer array of theoretical ion masses for all candidates flattened.</param>
/// <param name="candidatesIdx">An integer array that contains indices of where each candidate starts in candidatesValues.</param>
/// <param name="spectraValues">An integer array of peaks from experimental spectra flattened.</param>
/// <param name="spectraIdx">An integer array that contains indices of where each spectrum starts in spectraValues.</param>
/// <param name="cVLength">Length (int) of candidatesValues.</param>
/// <param name="cILength">Length (int) of candidatesIdx.</param>
/// <param name="sVLength">Length (int) of spectraValues.</param>
/// <param name="sILength">Length (int) of spectraIdx.</param>
/// <param name="n">How many of the best hits should be returned (int).</param>
/// <param name="tolerance">Tolerance for peak matching (float).</param>
/// <param name="normalize">If candidate vectors should be normalized to sum(elements) = 1 (bool).</param>
/// <param name="gaussianTol">If spectrum peaks should be modelled as normal distributions or not (bool).</param>
/// <param name="cores">Number of cores (int) used.</param>
/// <param name="verbose">Print info every (int) processed spectra.</param>
/// <returns>An integer array of length sILength * n containing the indexes of the top n candidates for each spectrum.</returns>
int* findTopCandidates(int* candidatesValues, int* candidatesIdx, 
                       int* spectraValues, int* spectraIdx, 
                       int cVLength, int cILength, 
                       int sVLength, int sILength,
                       int n, float tolerance,
                       bool normalize, bool gaussianTol,
                       int cores, int verbose) {

    if (n > cILength) {
        throw std::invalid_argument("Cannot return more hits than number of candidates!");
    }

    std::cout << "Running custom search version " << versionMajor << "." << versionMinor << "." << versionFix << std::endl;

    /// #### MATRIX CREATION ####
    // CSR nomenclature as used by cuSPARSE: https://docs.nvidia.com/cuda/cusparse/index.html#compressed-sparse-row-csr
    // the matrix in CSR format is given like this:
    int numberRows = cILength;
    const int numberColumns = ENCODING_SIZE;
    std::vector<int> rowOffsets;
    std::vector<int> columnIndices;
    std::vector<float> values;

    // create the row offsets vector
    rowOffsets.reserve(cILength + 1);
    for (int i = 0; i < cILength; ++i) {
        rowOffsets.push_back(candidatesIdx[i]);
    }
    rowOffsets.push_back(cILength);

    // create the column indices vector
    columnIndices.reserve(cVLength);
    for (int i = 0; i < cVLength; ++i) {
        columnIndices.push_back(candidatesValues[i]);
    }

    // create the values vector
    values.reserve(cVLength);
    for (int i = 0; i < cILength; ++i) {
        int startIter = candidatesIdx[i];
        int endIter = i + 1 == cILength ? cVLength : candidatesIdx[i + 1];
        int nrNonZero = endIter - startIter;
        float val = normalize ? 1.0 / (float) nrNonZero : 1.0;
        for (int j = startIter; j < endIter; ++j) {
            values.push_back(val);
        }
    }
    // #### END OF MATRIX CREATION ####

    auto* result = new int[sILength * n];
    float t = round(tolerance * MASS_MULTIPLIER);

    // for SPMV we have to iterate over every encoding vector
    for (int i = 0; i < sILength; ++i) {

        // #### VECTOR CREATION ####
        int startIter = spectraIdx[i];
        int endIter = i + 1 == sILength ? sVLength : spectraIdx[i + 1];
        std::vector<float> v(ENCODING_SIZE, 0.0f);
        for (int j = startIter; j < endIter; ++j) {
            auto currentPeak = spectraValues[j];
            auto minPeak = currentPeak - t > 0 ? currentPeak - t : 0;
            auto maxPeak = currentPeak + t < ENCODING_SIZE ? currentPeak + t : ENCODING_SIZE - 1;

            for (int k = minPeak; k <= maxPeak; ++k) {
                float newVal = gaussianTol ? normpdf((float) k, (float) currentPeak, (float) (t / 3.0)) : 1.0;
                v[k] = max(v[k], newVal);
            }
        }
        // #### END OF VECTOR CREATION ####

        // #### MULTIPLICATION ####
        // this is an examplary, very simple, non-optimized, single-threaded implementation of SPMV
        // for multi-threaded implementations the number of CPU cores used should be limited to a maximum of param "cores"
        std::vector<float> resultSPMV(numberRows, 0.0f);

        for (int row = 0; row < numberRows; ++row) {
            int rowBegin = rowOffsets[row];
            int rowEnd = rowOffsets[row + 1];
            float currentSum = 0.0f;
            for (int col = rowBegin; col < rowEnd; ++col) {
                int colIdx = columnIndices[col];
                float val = values[col];
                currentSum += (val * v[colIdx]);
            }
            resultSPMV[row] = currentSum;
        }
        // #### END OF MULTIPLICATION ####

        // #### GET MAX INDICES ####
        auto* idx = new int[numberRows];
        std::iota(idx, idx + numberRows, 0);
        std::sort(idx, idx + numberRows, [&](int i, int j) {return resultSPMV[i] > resultSPMV[j];});

        // assigning values to final result array
        for (int j = 0; j < n; ++j) {
            result[i * n + j] = idx[j];
        }

        delete[] idx;
        // #### END OF GET MAX INDICES ####

        if (verbose != 0 && (i + 1) % verbose == 0) {
            std::cout << "Searched " << i + 1 << " spectra in total..." << std::endl;
        }
    }

    return result;
}

#pragma region NotImplementedI

/// <summary>
/// A function that calculates the top n candidates for each spectrum. Maps to CPU_SVi32 (and variations).
/// </summary>
/// <param name="candidatesValues">An integer array of theoretical ion masses for all candidates flattened.</param>
/// <param name="candidatesIdx">An integer array that contains indices of where each candidate starts in candidatesValues.</param>
/// <param name="spectraValues">An integer array of peaks from experimental spectra flattened.</param>
/// <param name="spectraIdx">An integer array that contains indices of where each spectrum starts in spectraValues.</param>
/// <param name="cVLength">Length (int) of candidatesValues.</param>
/// <param name="cILength">Length (int) of candidatesIdx.</param>
/// <param name="sVLength">Length (int) of spectraValues.</param>
/// <param name="sILength">Length (int) of spectraIdx.</param>
/// <param name="n">How many of the best hits should be returned (int).</param>
/// <param name="tolerance">Tolerance for peak matching (float).</param>
/// <param name="normalize">If candidate vectors should be normalized to sum(elements) = 1 (bool).</param>
/// <param name="gaussianTol">If spectrum peaks should be modelled as normal distributions or not (bool).</param>
/// <param name="cores">Number of cores (int) used.</param>
/// <param name="verbose">Print info every (int) processed spectra.</param>
/// <returns>An integer array of length sILength * n containing the indexes of the top n candidates for each spectrum.</returns>
int* findTopCandidatesInt(int* candidatesValues, int* candidatesIdx,
                          int* spectraValues, int* spectraIdx,
                          int cVLength, int cILength,
                          int sVLength, int sILength,
                          int n, float tolerance,
                          bool normalize, bool gaussianTol,
                          int cores, int verbose) {

    if (n > cILength) {
        throw std::invalid_argument("Cannot return more hits than number of candidates!");
    }

    if (tolerance < 0.01f) {
        throw std::invalid_argument("Tolerance must not be smaller than 0.01 for i32 operations!");
    }

    throw NotImplementedException();

    auto* result = new int[sILength * n] {0};

    return result;
}

/// <summary>
/// A function that calculates the top n candidates for each spectrum. Maps to CPU_DVf32 (and variations). 
/// </summary>
/// <param name="candidatesValues">An integer array of theoretical ion masses for all candidates flattened.</param>
/// <param name="candidatesIdx">An integer array that contains indices of where each candidate starts in candidatesValues.</param>
/// <param name="spectraValues">An integer array of peaks from experimental spectra flattened.</param>
/// <param name="spectraIdx">An integer array that contains indices of where each spectrum starts in spectraValues.</param>
/// <param name="cVLength">Length (int) of candidatesValues.</param>
/// <param name="cILength">Length (int) of candidatesIdx.</param>
/// <param name="sVLength">Length (int) of spectraValues.</param>
/// <param name="sILength">Length (int) of spectraIdx.</param>
/// <param name="n">How many of the best hits should be returned (int).</param>
/// <param name="tolerance">Tolerance for peak matching (float).</param>
/// <param name="normalize">If candidate vectors should be normalized to sum(elements) = 1 (bool).</param>
/// <param name="gaussianTol">If spectrum peaks should be modelled as normal distributions or not (bool).</param>
/// <param name="cores">Number of cores (int) used.</param>
/// <param name="verbose">Print info every (int) processed spectra.</param>
/// <returns>An integer array of length sILength * n containing the indexes of the top n candidates for each spectrum.</returns>
int* findTopCandidates2(int* candidatesValues, int* candidatesIdx,
                        int* spectraValues, int* spectraIdx,
                        int cVLength, int cILength,
                        int sVLength, int sILength,
                        int n, float tolerance,
                        bool normalize, bool gaussianTol,
                        int cores, int verbose) {

    if (n > cILength) {
        throw std::invalid_argument("Cannot return more hits than number of candidates!");
    }

    throw NotImplementedException();

    auto* result = new int[sILength * n] {0};

    return result;
}

/// <summary>
/// A function that calculates the top n candidates for each spectrum. Maps to CPU_DVi32 (and variations).
/// </summary>
/// <param name="candidatesValues">An integer array of theoretical ion masses for all candidates flattened.</param>
/// <param name="candidatesIdx">An integer array that contains indices of where each candidate starts in candidatesValues.</param>
/// <param name="spectraValues">An integer array of peaks from experimental spectra flattened.</param>
/// <param name="spectraIdx">An integer array that contains indices of where each spectrum starts in spectraValues.</param>
/// <param name="cVLength">Length (int) of candidatesValues.</param>
/// <param name="cILength">Length (int) of candidatesIdx.</param>
/// <param name="sVLength">Length (int) of spectraValues.</param>
/// <param name="sILength">Length (int) of spectraIdx.</param>
/// <param name="n">How many of the best hits should be returned (int).</param>
/// <param name="tolerance">Tolerance for peak matching (float >= 0.01).</param>
/// <param name="normalize">If candidate vectors should be normalized to sum(elements) = 1 (bool).</param>
/// <param name="gaussianTol">If spectrum peaks should be modelled as normal distributions or not (bool).</param>
/// <param name="cores">Number of cores (int) used.</param>
/// <param name="verbose">Print info every (int) processed spectra.</param>
/// <returns>An integer array of length sILength * n containing the indexes of the top n candidates for each spectrum.</returns>
int* findTopCandidates2Int(int* candidatesValues, int* candidatesIdx,
                           int* spectraValues, int* spectraIdx,
                           int cVLength, int cILength,
                           int sVLength, int sILength,
                           int n, float tolerance,
                           bool normalize, bool gaussianTol,
                           int cores, int verbose) {

    if (n > cILength) {
        throw std::invalid_argument("Cannot return more hits than number of candidates!");
    }

    if (tolerance < 0.01f) {
        throw std::invalid_argument("Tolerance must not be smaller than 0.01 for i32 operations!");
    }

    throw NotImplementedException();

    auto* result = new int[sILength * n] {0};

    return result;
}

#pragma endregion NotImplementedI

/// <summary>
/// A function that calculates the top n candidates for each spectrum. Maps to CPU_SMf32 (and variations).
/// </summary>
/// <param name="candidatesValues">An integer array of theoretical ion masses for all candidates flattened.</param>
/// <param name="candidatesIdx">An integer array that contains indices of where each candidate starts in candidatesValues.</param>
/// <param name="spectraValues">An integer array of peaks from experimental spectra flattened.</param>
/// <param name="spectraIdx">An integer array that contains indices of where each spectrum starts in spectraValues.</param>
/// <param name="cVLength">Length (int) of candidatesValues.</param>
/// <param name="cILength">Length (int) of candidatesIdx.</param>
/// <param name="sVLength">Length (int) of spectraValues.</param>
/// <param name="sILength">Length (int) of spectraIdx.</param>
/// <param name="n">How many of the best hits should be returned (int).</param>
/// <param name="tolerance">Tolerance for peak matching (float).</param>
/// <param name="normalize">If candidate vectors should be normalized to sum(elements) = 1 (bool).</param>
/// <param name="gaussianTol">If spectrum peaks should be modelled as normal distributions or not (bool).</param>
/// <param name="batchSize">How many spectra (int) should be searched at once.</param>
/// <param name="cores">Number of cores (int) used.</param>
/// <param name="verbose">Print info every (int) processed spectra.</param>
/// <returns>An integer array of length sILength * n containing the indexes of the top n candidates for each spectrum.</returns>
int* findTopCandidatesBatched(int* candidatesValues, int* candidatesIdx,
                              int* spectraValues, int* spectraIdx,
                              int cVLength, int cILength,
                              int sVLength, int sILength,
                              int n, float tolerance,
                              bool normalize, bool gaussianTol,
                              int batchSize,
                              int cores, int verbose) {

    if (n > cILength) {
        throw std::invalid_argument("Cannot return more hits than number of candidates!");
    }

    // #### SPMM ####
    // consider a SPMM: A * B = C
    // matrix A is a sparse matrix as given in "findTopCandidates"
    // matrix B and C are sparse or dense, this can be chosen depending on the needs
    // for SPMM we have the param "batchSize" that determines the size of intermediate matrices B for multiplication
    // this is due to usually not being able to process all values of params "spectraValues" and "spectraIdx" at
    // once
    // for more information, please refer to how this is handled for dense matrices here:
    // https://github.com/hgb-bin-proteomics/CandidateVectorSearch/blob/57fc3b42a32df0405e3fbf08c908ec780510417a/VectorSearch/dllmain.cpp#L852
    // and for sparse matrices here:
    // https://github.com/hgb-bin-proteomics/CandidateVectorSearch/blob/57fc3b42a32df0405e3fbf08c908ec780510417a/VectorSearch/dllmain.cpp#L570

    throw NotImplementedException();

    auto* result = new int[sILength * n] {0};

    return result;
}

#pragma region NotImplementedII

/// <summary>
/// A function that calculates the top n candidates for each spectrum. Maps to CPU_SMi32 (and variations).
/// </summary>
/// <param name="candidatesValues">An integer array of theoretical ion masses for all candidates flattened.</param>
/// <param name="candidatesIdx">An integer array that contains indices of where each candidate starts in candidatesValues.</param>
/// <param name="spectraValues">An integer array of peaks from experimental spectra flattened.</param>
/// <param name="spectraIdx">An integer array that contains indices of where each spectrum starts in spectraValues.</param>
/// <param name="cVLength">Length (int) of candidatesValues.</param>
/// <param name="cILength">Length (int) of candidatesIdx.</param>
/// <param name="sVLength">Length (int) of spectraValues.</param>
/// <param name="sILength">Length (int) of spectraIdx.</param>
/// <param name="n">How many of the best hits should be returned (int).</param>
/// <param name="tolerance">Tolerance for peak matching (float).</param>
/// <param name="normalize">If candidate vectors should be normalized to sum(elements) = 1 (bool).</param>
/// <param name="gaussianTol">If spectrum peaks should be modelled as normal distributions or not (bool).</param>
/// <param name="batchSize">How many spectra (int) should be searched at once.</param>
/// <param name="cores">Number of cores (int) used.</param>
/// <param name="verbose">Print info every (int) processed spectra.</param>
/// <returns>An integer array of length sILength * n containing the indexes of the top n candidates for each spectrum.</returns>
int* findTopCandidatesBatchedInt(int* candidatesValues, int* candidatesIdx,
                                 int* spectraValues, int* spectraIdx,
                                 int cVLength, int cILength,
                                 int sVLength, int sILength,
                                 int n, float tolerance,
                                 bool normalize, bool gaussianTol,
                                 int batchSize,
                                 int cores, int verbose) {

    if (n > cILength) {
        throw std::invalid_argument("Cannot return more hits than number of candidates!");
    }

    if (tolerance < 0.01f) {
        throw std::invalid_argument("Tolerance must not be smaller than 0.01 for i32 operations!");
    }

    throw NotImplementedException();

    auto* result = new int[sILength * n] {0};

    return result;
}

/// <summary>
/// A function that calculates the top n candidates for each spectrum. Maps to CPU_DMf32 (and variations).
/// </summary>
/// <param name="candidatesValues">An integer array of theoretical ion masses for all candidates flattened.</param>
/// <param name="candidatesIdx">An integer array that contains indices of where each candidate starts in candidatesValues.</param>
/// <param name="spectraValues">An integer array of peaks from experimental spectra flattened.</param>
/// <param name="spectraIdx">An integer array that contains indices of where each spectrum starts in spectraValues.</param>
/// <param name="cVLength">Length (int) of candidatesValues.</param>
/// <param name="cILength">Length (int) of candidatesIdx.</param>
/// <param name="sVLength">Length (int) of spectraValues.</param>
/// <param name="sILength">Length (int) of spectraIdx.</param>
/// <param name="n">How many of the best hits should be returned (int).</param>
/// <param name="tolerance">Tolerance for peak matching (float).</param>
/// <param name="normalize">If candidate vectors should be normalized to sum(elements) = 1 (bool).</param>
/// <param name="gaussianTol">If spectrum peaks should be modelled as normal distributions or not (bool).</param>
/// <param name="batchSize">How many spectra (int) should be searched at once.</param>
/// <param name="cores">Number of cores (int) used.</param>
/// <param name="verbose">Print info every (int) processed spectra.</param>
/// <returns>An integer array of length sILength * n containing the indexes of the top n candidates for each spectrum.</returns>
int* findTopCandidatesBatched2(int* candidatesValues, int* candidatesIdx,
                               int* spectraValues, int* spectraIdx,
                               int cVLength, int cILength,
                               int sVLength, int sILength,
                               int n, float tolerance,
                               bool normalize, bool gaussianTol,
                               int batchSize,
                               int cores, int verbose) {

    if (n > cILength) {
        throw std::invalid_argument("Cannot return more hits than number of candidates!");
    }

    throw NotImplementedException();

    auto* result = new int[sILength * n] {0};
  
    return result;
}

/// <summary>
/// A function that calculates the top n candidates for each spectrum. Maps to CPU_DMi32 (and variations).
/// </summary>
/// <param name="candidatesValues">An integer array of theoretical ion masses for all candidates flattened.</param>
/// <param name="candidatesIdx">An integer array that contains indices of where each candidate starts in candidatesValues.</param>
/// <param name="spectraValues">An integer array of peaks from experimental spectra flattened.</param>
/// <param name="spectraIdx">An integer array that contains indices of where each spectrum starts in spectraValues.</param>
/// <param name="cVLength">Length (int) of candidatesValues.</param>
/// <param name="cILength">Length (int) of candidatesIdx.</param>
/// <param name="sVLength">Length (int) of spectraValues.</param>
/// <param name="sILength">Length (int) of spectraIdx.</param>
/// <param name="n">How many of the best hits should be returned (int).</param>
/// <param name="tolerance">Tolerance for peak matching (float >= 0.01).</param>
/// <param name="normalize">If candidate vectors should be normalized to sum(elements) = 1 (bool).</param>
/// <param name="gaussianTol">If spectrum peaks should be modelled as normal distributions or not (bool).</param>
/// <param name="batchSize">How many spectra (int) should be searched at once.</param>
/// <param name="cores">Number of cores (int) used.</param>
/// <param name="verbose">Print info every (int) processed spectra.</param>
/// <returns>An integer array of length sILength * n containing the indexes of the top n candidates for each spectrum.</returns>
int* findTopCandidatesBatched2Int(int* candidatesValues, int* candidatesIdx,
                                  int* spectraValues, int* spectraIdx,
                                  int cVLength, int cILength,
                                  int sVLength, int sILength,
                                  int n, float tolerance,
                                  bool normalize, bool gaussianTol,
                                  int batchSize,
                                  int cores, int verbose) {

    if (n > cILength) {
        throw std::invalid_argument("Cannot return more hits than number of candidates!");
    }

    if (tolerance < 0.01f) {
        throw std::invalid_argument("Tolerance must not be smaller than 0.01 for i32 operations!");
    }

    throw NotImplementedException();

    auto* result = new int[sILength * n] {0};

    return result;
}

#pragma endregion NotImplementedII

/// <summary>
/// Free memory after result has been marshalled.
/// </summary>
/// <param name="result">The result array.</param>
/// <returns>0</returns>
int releaseMemory(int* result) {

    delete[] result;
    return 0;
}

/// <summary>
/// Returns the square for a given value x.
/// </summary>
/// <param name="x">The value to be squared.</param>
/// <returns>Square of x.</returns>
float squared(float x) {
    return x * x;
}

/// <summary>
/// Returns the PDF for a given x for the normal distribution given by mu and sigma.
/// </summary>
/// <param name="x">The value for which the PDF should be calculated.</param>
/// <param name="mu">The mu of the normal distribution.</param>
/// <param name="sigma">The sigma of the normal distribution.</param>
/// <returns>The PDF at x for the normal distribution given by mu and sigma. If sigma = 0 it returns 1.</returns>
float normpdf(float x, float mu, float sigma) {
    if (sigma == 0.0) {
        return 1.0;
    }
    return (ONE_OVER_SQRT_PI / sigma) * exp(-0.5 * squared((x - mu) / sigma));
}

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
    switch (ul_reason_for_call)
    {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
        break;
    }
    return TRUE;
}

