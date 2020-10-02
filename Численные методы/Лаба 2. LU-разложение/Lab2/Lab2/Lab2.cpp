#include <iostream>
#include <algorithm>

using namespace std;

const double minValue = 0.00001;

void createMatrix(double*** pMatrix, int N) {
    auto& matrix = *pMatrix;
    matrix = new double* [N];
    for (int i = 0; i < N; i++) {
        matrix[i] = new double[N];
    }
}

template<int N>
void fillMatrix(double** matrix, double values[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = values[i][j];
        }
    }
}

void fillMatrixAsEmpty(double** matrix, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = 0.0;
        }
    }
}

void printMatrix(double** matrix, int N) {
    for (int i = 0; i < N; i++) {
        printf("| ");
        for (int j = 0; j < N; j++) {
            printf("%.2f ", matrix[i][j]);
        }
        printf("|\n");
    }
    printf("\n");
}

void matrixMul(double** A, double** B, double** C, int N)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                C[i][j] += A[i][k] * B[k][j];
}

int defineRowIdxWithMainValue(double** matrix, int k, int N) {
    int m = k;
    double maxValue = matrix[m][k];
    for (int i = k + 1; i < N; i++) {
        if (matrix[i][k] > maxValue) {
            m = i;
            maxValue = matrix[m][k];
        }
    }
    return m;
}

void LUdecomposition(double** L, double** U, int& rank, int N) {
    for (int k = 0; k < N; k++) {
        auto rowIdx = defineRowIdxWithMainValue(U, k, N);
        swap(U[k], U[rowIdx]);

        double mainValue = U[k][k];
        if (mainValue < minValue) {
            rank = k;
            return;
        }

        for (int j = k; j < N; j++) {
            U[k][j] /= mainValue;
            //printMatrix(U, N);
        }

        for (int i = k + 1; i < N; i++) {
            auto val = U[i][k];
            L[i][k] = val / mainValue;
            for (int j = k; j < N; j++) {
                U[i][j] -= val * U[k][j];
                //printMatrix(U, N);
            }
        }
        printMatrix(U, N);
    }
}

int main()
{
    /*const int N = 4;
    double values[N][N] = {
        { 2.0, 3.0, 4.0, 5.0 },
        { 2.0, 3.0, 4.0, 6.0 },
        { 2.0, 3.0, 4.0, 7.0 },
        { 2.0, 3.0, 4.0, 8.0 }
    };*/

    const int N = 3;
    double values[N][N] = {
        { 10, -7, 0 },
        {-3, 6, 2},
        {5, -1, 5}
    };

    double** A = nullptr;
    createMatrix(&A, N);
    fillMatrix<N>(A, values);

    double** U = nullptr;
    createMatrix(&U, N);
    fillMatrix<N>(U, values);

    double** L = nullptr;
    createMatrix(&L, N);
    fillMatrixAsEmpty(L, N);

    printf("Matrix A:\n");
    printMatrix(A, N);
    int rank = N;
    LUdecomposition(L, U, rank, N);

    printf("Matrix L:\n");
    printMatrix(L, N);
    printf("Matrix U:\n");
    printMatrix(U, N);
    printf("\n\nRank = %i\n", rank);


    double** C = nullptr;
    createMatrix(&C, N);
    fillMatrixAsEmpty(C, N);
    matrixMul(L, U, C, N);
    printMatrix(C, N);
}