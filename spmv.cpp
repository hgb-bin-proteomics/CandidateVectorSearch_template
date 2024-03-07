#include <vector>
#include <iostream>

int main() {

    /**
     * PROBLEM DEFINITION
     *
     * SPMV = A * V
     * A: sparse matrix
     * V: dense vector
     * SPMV: dense vector
     *
     * MATRIX A
     * [[1, 0, 2, 3],
     *  [0, 0, 0, 0],
     *  [0, 4, 0, 0],
     *  [5, 0, 6, 7],
     *  [0, 8, 0, 9]]
     *
     * VECTOR V
     * [1, 2, 3, 4]
    **/

    // MATRIX A
    std::vector<int> rowOffsets{0, 3, 3, 4, 7, 9};
    std::vector<int> columnIndices{0, 2, 3, 1, 0, 2, 3, 1, 3};
    std::vector<float> values{1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f};
    int numberRows = rowOffsets.size() - 1;

    // VECTOR V
    std::vector<float> v{1.0f, 2.0f, 3.0f, 4.0f};

    // RESULT VECTOR
    std::vector<float> resultSPMV(numberRows, 0.0f);

    // MULTIPLICATION
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

    // PRINT RESULT
    for (int i = 0; i < numberRows; ++i) {
        std::cout << resultSPMV[i] << std::endl;
    }

    return 0;
}