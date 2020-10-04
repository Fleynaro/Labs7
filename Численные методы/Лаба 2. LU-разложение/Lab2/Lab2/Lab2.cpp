#include <iostream>
#include <algorithm>
#define M_PI            3.14159265358979323846

using namespace std;

const double minValue = 0.00001;

//создаем в памяти матрицу
void createMatrix(double*** pMatrix, int N) {
    auto& matrix = *pMatrix;
    matrix = new double* [N];
    for (int i = 0; i < N; i++) {
        matrix[i] = new double[N];
    }
}

//заполняем матрицу значениями
template<int N>
void fillMatrix(double** matrix, double values[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = values[i][j];
        }
    }
}

//заполняем матрицу нулями
void fillMatrixAsEmpty(double** matrix, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = 0.0;
        }
    }
}

//копируем значения одной матрицы в другую
void copyMatrixToMatrix(double** srcMatrix, double** dstMatrix, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            dstMatrix[i][j] = srcMatrix[i][j];
        }
    }
}

//распечатка вектора
void printVector(double* vector, int N) {
    printf("| ");
    for (int i = 0; i < N; i++) {
        printf("%.2f ", vector[i]);
    }
    printf("|\n");
}

//распечатка матрицы
void printMatrix(double** matrix, int N) {
    for (int i = 0; i < N; i++) {
        printVector(matrix[i], N);
    }
    printf("\n");
}

//умножение матрицы
void matrixMul(double** A, double** B, double** C, int N)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                C[i][j] += A[i][k] * B[k][j];
}

//определить строку с главным элементом (который максимален в текущем столбце)
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

//LU разложение
void LUdecomposition(double** L, double** U, int* P, int& rank, double& sign, int N) {
    for (int k = 0; k < N; k++) {
        auto rowIdx = defineRowIdxWithMainValue(U, k, N);
        if (k != rowIdx) {
            //Смена строк
            swap(U[k], U[rowIdx]);
            P[k] = rowIdx;
            P[rowIdx] = k;
            sign *= -1.0;
        }
        else {
            P[k] = k;
        }

        //главный элемент
        double mainValue = U[k][k];
        if (mainValue < minValue) {
            //Определяем ранг матрицы
            rank = k;
            return;
        }

        //заполняем матрицу L
        for (int i = k; i < N; i++) {
            L[i][k] = U[i][k] / mainValue;
            //printMatrix(L, N);
        }

        //заполняем матрицу U
        for (int i = k + 1; i < N; i++) {
            for (int j = k; j < N; j++) {
                U[i][j] = U[i][j] - L[i][k] * U[k][j];
                //printMatrix(U, N);
            }
        }
        //printMatrix(U, N);
    }
}


//Решение уравнения Ly = Pb
void SolveLy(double** triangleMatrix, double* X, double* B, int N) {
    for (int i = 0; i < N; i++) {
        X[i] = B[i];
        for (int j = 0; j < i; j++) {
            X[i] -= X[j] * triangleMatrix[i][j];
        }
    }
}

//Решение уравнения Ux = y
void SolveUx(double** triangleMatrix, double* X, double* B, int N) {
    for (int i = N - 1; i >= 0; i--) {
        X[i] = B[i] / triangleMatrix[i][i];
        for (int j = N - 1; j > i; j--) {
            X[i] -= X[j] * triangleMatrix[i][j] / triangleMatrix[i][i];
        }
    }
}

//Решение уравнения Ax = b, то есть LUx = Pb
void SolveSOLE(double** L, double** U, double* X, int* P, double* B, int N) {
    //Ax = b (A = PLU)
    //LUx = Pb (Ux = y)
    //Ly = Pb
    double* vectorY = new double[N];
    //"умножаем" вектор b на матрицу перестановок
    double* vectorPB = new double[N];
    for (int i = 0; i < N; i++) {
        vectorPB[i] = B[P[i]];
    }
    //printVector(PB, N);
    SolveLy(L, vectorY, vectorPB, N);
    //printVector(Y, N);
    SolveUx(U, X, vectorY, N);
}

void SolveBackwardMatrix(double** L, double** U, double** X, int* P, int N) {
    //LUX = PE
    double* vectorX = new double[N];
    double* vectorE = new double[N];
    for (int t = 0; t < N; t++) {
        vectorE[t] = 0.0;
    }
    for (int i = 0; i < N; i++) {
        //формируем вектор Ei
        if (i != 0) {
            vectorE[i - 1] = 0.0;
        }
        vectorE[i] = 1.0;

        //получаем вектор-столбец X
        SolveSOLE(L, U, vectorX, P, vectorE, N);
        //записываем его в матрицу X
        for (int t = 0; t < N; t++) {
            X[t][i] = vectorX[t];
        }
    }
}

//найти определитель
double computeDet(double** U, int N, double sign) {
    double det = sign;
    for (int t = 0; t < N; t++) {
        det *= U[t][t];
    }
    return det;
}

//транспонировать матрицу
void transpose(double** matrix, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            swap(matrix[j][i], matrix[i][j]);
        }
    }
}

//вычислить кубическую норму
double computeCubNorm(double** matrix, int N) {
    double result = 0.0;
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < N; j++) {
            sum += abs(matrix[i][j]);
        }
        result = max(result, sum);
    }
    return result;
}
//////////////////////////////////////////////////////////////////////////////////////////////////

//Нахождение максимального элемента в матрице, находящегося не на диагонали
void searchMaxElemMatrix(double** matrix, const int N, int& imax, int& jmax) {
    double max = 0.0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i != j && abs(matrix[i][j]) > max)
            {
                max = abs(matrix[i][j]);
                imax = i;
                jmax = j;
            }
        }
    }
}

//tan(2*alpha) = 2*a[i][j]/(a[i][i]-a[j][j])
double getAlpha(double** matrix, int imax, int jmax) {
    double alpha;
    if (matrix[imax][imax] - matrix[jmax][jmax] == 0)
    {
        alpha = M_PI / 4;
    }
    else
    {
        alpha = atan(2 * matrix[imax][jmax] / (matrix[imax][imax] - matrix[jmax][jmax])) / 2;
    }
    return alpha;
}

//Является ли матрица диагональной
bool isMatrixDiagonal(double** matrix, const int N) {
    double kvSum = 0.0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i != j)
            {
                kvSum += matrix[i][j] * matrix[i][j];
            }
        }
    }
    return kvSum < 2.0;
}

//Получение диагональной матрицы
double** getNewDiagonalMatrixByRotation(double** matrix, const int N) {
    double** rotatedMatrix = nullptr;
    createMatrix(&rotatedMatrix, N);
    copyMatrixToMatrix(matrix, rotatedMatrix, N);

    //вспомогательная матрица
    double** B = nullptr;
    createMatrix(&B, N);

    while (!isMatrixDiagonal(rotatedMatrix, N)) {
        printMatrix(rotatedMatrix, N);
        int imax, jmax;
        searchMaxElemMatrix(rotatedMatrix, N, imax, jmax);
        double alpha = getAlpha(rotatedMatrix, imax, jmax);
        double c = cos(alpha);
        double s = sin(alpha);

        //результат умножения матрицы A в k-ом состоянии на матрицу вращения справа
        copyMatrixToMatrix(rotatedMatrix, B, N);
        for (int m = 0; m < N; m++) {
            B[m][imax] = c * rotatedMatrix[m][imax] + s * rotatedMatrix[m][jmax];
            B[m][jmax] = -s * rotatedMatrix[m][imax] + c * rotatedMatrix[m][jmax];
        }

        //результат умножения матрицы B на матрицу вращения слева
        for (int m = 0; m < N; m++) {
            rotatedMatrix[imax][m] = c * B[imax][m] + s * B[jmax][m];
            rotatedMatrix[jmax][m] = -s * B[imax][m] + c * B[jmax][m];
        }
    }

    return rotatedMatrix;
}

//Получение максимального собственного значения
double getMaxEigenvalue(double** matrix, const int N) {
    double max = 0;
    for (int i = 0; i < N; i++)
    {
        if (abs(matrix[i][i]) > max)
            max = abs(matrix[i][i]);
    }
    return max;
}

//Вычисление евклидовой нормы матрицы
double computEuclidNorm(double** A, double** trA, const int N) {
    double** newMatrix = nullptr;
    createMatrix(&newMatrix, N);
    fillMatrixAsEmpty(newMatrix, N);
    matrixMul(trA, A, newMatrix, N);
    printMatrix(newMatrix, N);
    newMatrix = getNewDiagonalMatrixByRotation(newMatrix, N);
    printMatrix(newMatrix, N);
    double eigenvalue = getMaxEigenvalue(newMatrix, N);
    return sqrt(eigenvalue);
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
    double matrixValues[N][N] = {
        { 10, -7, 0 },
        {-3, 6, 2},
        {5, -1, 5}
    };
    double vectorB[N] = {
        1, 2, 3
    };

    //Матрица A
    double** A = nullptr;
    createMatrix(&A, N);
    fillMatrix<N>(A, matrixValues);

    //Матрица U
    double** U = nullptr;
    createMatrix(&U, N);
    copyMatrixToMatrix(A, U, N);

    //Матрица L
    double** L = nullptr;
    createMatrix(&L, N);
    fillMatrixAsEmpty(L, N);

    printf("Matrix A:\n");
    printMatrix(A, N);
    printf("Vector B:\n");
    printVector(vectorB, N);

    //LU разложение
    int rank = N; //ранг матрицы
    int P[N]; //"матрица" перестановок (на самом деле подстановка)
    double sign = 1.0;
    LUdecomposition(L, U, P, rank, sign, N);

    printf("Matrix L:\n");
    printMatrix(L, N);
    printf("Matrix U:\n");
    printMatrix(U, N);
    printf("\n\nRank = %i\n", rank);


    double** C = nullptr;
    createMatrix(&C, N);
    fillMatrixAsEmpty(C, N);
    matrixMul(L, U, C, N);
    printf("Matrix LU:\n");
    printMatrix(C, N);


    //решаем СЛАУ
    double vectorX[N];
    SolveSOLE(L, U, vectorX, P, vectorB, N);
    printf("Vector X:\n");
    printVector(vectorX, N);


    //находим обратную матрицу
    double** backwardMatrix = nullptr;
    createMatrix(&backwardMatrix, N);
    SolveBackwardMatrix(L, U, backwardMatrix, P, N);
    printf("Backward matrix:\n");
    printMatrix(backwardMatrix, N);

    //найти определитель
    auto det = computeDet(U, N, sign);
    printf("det = %f\n\n", det);

    //транспонируем матрицу
    double** trA = nullptr;
    createMatrix(&trA, N);
    copyMatrixToMatrix(A, trA, N);
    transpose(trA, N);
    //транспонируем обратную матрицу
    double** trBackwardMatrix = nullptr;
    createMatrix(&trBackwardMatrix, N);
    copyMatrixToMatrix(backwardMatrix, trBackwardMatrix, N);
    transpose(trBackwardMatrix, N);


   /* double sMatrixValues[3][3] = {
        { 10, 4, 5 },
        {4, 20, 1},
        {5, 1, 30}
    };
    double** sMatrix = nullptr;
    createMatrix(&sMatrix, N);
    fillMatrix<3>(sMatrix, sMatrixValues);
    printMatrix(sMatrix, N);
    auto diagonalMatrix = getNewDiagonalMatrixByRotation(sMatrix, N);
    printMatrix(diagonalMatrix, N);*/

    //вычисление норм матриц
    auto cubNorm1 = computeCubNorm(A, N);
    auto cubNorm2 = computeCubNorm(backwardMatrix, N);
    auto octNorm1 = computeCubNorm(trA, N);
    auto octNorm2 = computeCubNorm(trBackwardMatrix, N);
    ////////////////////////////////////////////////////////
    double euclidNorm1 = computEuclidNorm(A, trA, N);
    double euclidNorm2 = computEuclidNorm(backwardMatrix, trBackwardMatrix, N);

    /////////////////////////////////////////////////////////

    //вычисление числа обусловленности
    auto cubCond = cubNorm1 * cubNorm2;
    auto octCond = octNorm1 * octNorm2;
    double euclidCond = euclidNorm1 * euclidNorm2;

    printf("cubCond = %f, octCond = %f, euclidCond = %f\n", cubCond, octCond, euclidCond);
    return 0;
}