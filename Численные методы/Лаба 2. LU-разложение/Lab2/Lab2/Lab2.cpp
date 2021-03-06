﻿#include <iostream>
#include <algorithm>
#define M_PI            3.14159265358979323846

using namespace std;
const double minValue = 0.00001;

//распечатка вектора
void printVector(double* vector, int N, bool expFormat = false) {
    printf(" ");
    for (int i = 0; i < N; i++) {
        if (expFormat) {
            printf("%.10e ", vector[i]);
        }
        else {
            printf("%.7f ", vector[i]);
        }
    }
    printf("\n");
}

//распечатка матрицы
void printMatrix(double** matrix, int N, bool expFormat = false) {
    for (int i = 0; i < N; i++) {
        printVector(matrix[i], N, expFormat);
    }
    printf("\n");
}

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

//заполняем матрицу как единичную
void fillMatrixAsE(double** matrix, int N) {
    fillMatrixAsEmpty(matrix, N);
    for (int i = 0; i < N; i++) {
        matrix[i][i] = 1.0;
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

//умножение матрицы на вектор
void matrixMulVec(double** A, double* B, double* C, int N)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            C[i] += A[i][j] * B[j];
}

//умножение матрицы на матрицу
void matrixMul(double** A, double** B, double** C, int N)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                C[i][j] += A[i][k] * B[k][j];
}

//вычитание матрицы
void matrixSub(double** A, double** B, double** C, int N)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            C[i][j] = A[i][j] - B[i][j];
}

//получение матрицы PA
double** matrixPA(double** A, int* P, int N)
{
    double** PA = new double* [N];
    for (int i = 0; i < N; i++) {
        PA[i] = A[P[i]];
    }
    return PA;
}

//вычитание векторов
void vectorSub(double* A, double* B, double* C, int N)
{
    for (int i = 0; i < N; i++)
        C[i] = A[i] - B[i];
}

//определить строку с главным элементом (который максимален в текущем столбце)
int defineRowIdxWithMainValue(double** matrix, int k, int N) {
    int m = k;
    double maxValue = 0.0;
    for (int i = k; i < N; i++) {
        if (abs(matrix[i][k]) > maxValue) {
            m = i;
            maxValue = abs(matrix[m][k]);
        }
    }
    return m;
}

//LU разложение
void LUdecomposition(double** L, double** U, int* P, int& rank, double& sign, int N) {
    //иницилизируем подстановку P
    for (int i = 0; i < N; i++) {
        P[i] = i;
    }

    printf("U = A:\n");
    printMatrix(U, N);
    printf("L:\n");
    printMatrix(L, N);

    for (int k = 0; k < N; k++) {
        auto rowIdx = defineRowIdxWithMainValue(U, k, N);
        if (k != rowIdx) {
            //Смена строк
            swap(U[k], U[rowIdx]);
            swap(L[k], L[rowIdx]);
            swap(P[k], P[rowIdx]);
            sign *= -1.0;
        }

        //главный элемент
        double mainValue = U[k][k];
        if (abs(mainValue) < minValue) {
            //Определяем ранг матрицы
            rank = k;
            return;
        }

        //заполняем матрицу L
        for (int i = k; i < N; i++) {
            L[i][k] = U[i][k];
            //printMatrix(L, N);
        }
        
        for (int j = k; j < N; j++) {
            U[k][j] /= mainValue;
        }
        
        //заполняем матрицу U
        for (int i = k + 1; i < N; i++) {
            for (int j = k; j < N; j++) {
                U[i][j] = U[i][j] - L[i][k] * U[k][j];
                //printMatrix(U, N);
            }
        }
        
        printf("\nk = %i\nm = %i\nU[m][k] = %.7f\nU:\n", k, rowIdx, mainValue);
        printMatrix(U, N);
        printf("L:\n");
        printMatrix(L, N);
    }
}


//Решение уравнения Ly = Pb
void SolveLy(double** triangleMatrix, double* X, double* B, int N) {
    for (int i = 0; i < N; i++) {
        X[i] = B[i] / triangleMatrix[i][i];
        for (int j = 0; j < i; j++) {
            X[i] -= X[j] * triangleMatrix[i][j] / triangleMatrix[i][i];
        }
    }
}

//Решение уравнения Ux = y
void SolveUx(double** triangleMatrix, double* X, double* B, int N) {
    for (int i = N - 1; i >= 0; i--) {
        X[i] = B[i];
        for (int j = N - 1; j > i; j--) {
            X[i] -= X[j] * triangleMatrix[i][j];
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
double computeDet(double** L, int N, double sign) {
    double det = sign;
    for (int t = 0; t < N; t++) {
        det *= L[t][t];
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
    return kvSum < minValue * minValue;
}

//Получение диагональной матрицы
double** getNewDiagonalMatrixByRotation(double** matrix, const int N) {
    double* vectorI = new double[N];
    double* vectorJ = new double[N];
    
    double** rotatedMatrix = nullptr;
    createMatrix(&rotatedMatrix, N);
    copyMatrixToMatrix(matrix, rotatedMatrix, N);

    //вспомогательная матрица
    double** B = nullptr;
    createMatrix(&B, N);

    while (!isMatrixDiagonal(rotatedMatrix, N)) {
        //printMatrix(rotatedMatrix, N);
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
            vectorI[m] = c * B[imax][m] + s * B[jmax][m];
            vectorJ[m] = -s * B[imax][m] + c * B[jmax][m];
        }

        swap(B[imax], vectorI);
        swap(B[jmax], vectorJ);
        copyMatrixToMatrix(B, rotatedMatrix, N);
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
    //printMatrix(newMatrix, N);
    newMatrix = getNewDiagonalMatrixByRotation(newMatrix, N);
    //printMatrix(newMatrix, N);
    double eigenvalue = getMaxEigenvalue(newMatrix, N);
    return sqrt(eigenvalue);
}

int main()
{
    const int N = 4;
    double matrixValues[N][N] = {
        { 0.1, -8.3, 7.1, 5.5},
        { 1.2, 5.2, -9.1, -0.2},
        {-7.9, 9.6, 0.9, -1.2},
        {3.8, -4.7, -0.2, 7.9}
    };
    double vectorB[N] = {
        26.8, -16.5, 9.2, 25.4
    };
    double vectorX_[N] = {
        1.0, 2.0, 3.0, 4.0
    };

    //единичная матрица
    double** E = nullptr;
    createMatrix(&E, N);
    fillMatrixAsE(E, N);

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

    printf("1) Input\A:\n");
    printMatrix(A, N);
    printf("B:\n");
    printVector(vectorB, N);

    //LU разложение
    int rank = N; //ранг матрицы
    int P[N]; //"матрица" перестановок (на самом деле подстановка)
    double sign = 1.0;

    printf("\n2) LU decomposition\n", rank);
    LUdecomposition(L, U, P, rank, sign, N);

    printf("\nRank = %i\n", rank);
    if (rank != N)
        return 0;

    auto PA = matrixPA(A, P, N);
    //LU матрица
    double** LU = nullptr;
    createMatrix(&LU, N);
    fillMatrixAsEmpty(LU, N);
    matrixMul(L, U, LU, N);
    if (false) { //show PA/LU
        printf("Matrix PA:\n");
        printMatrix(PA, N);
        printf("Matrix LU (check to see LU = PA):\n");
        printMatrix(LU, N);
    }
    
    //LU - PA
    printf("\nLU-PA:\n");
    double** LUminusPA = nullptr;
    createMatrix(&LUminusPA, N);
    matrixSub(LU, PA, LUminusPA, N);
    printMatrix(LUminusPA, N, true);

    //решаем СЛАУ
    double vectorX[N];
    SolveSOLE(L, U, vectorX, P, vectorB, N);
    printf("\n3) SOLE\nVector X:\n");
    printVector(vectorX, N);


    //вектор невязки
    double vectorX2[N];
    double vectorR[N];
    for (int i = 0; i < N; i++) vectorX2[i] = 0;
    matrixMulVec(A, vectorX, vectorX2, N);
    vectorSub(vectorX2, vectorB, vectorR, N);
    printf("\nAx - b = \n");
    printVector(vectorR, N, true);

    //погрешность решения
    double vectorDelta[N];
    vectorSub(vectorX_, vectorX, vectorDelta, N);
    printf("X* - X = \n");
    printVector(vectorDelta, N, true);


    //находим обратную матрицу
    double** backwardMatrix = nullptr;
    createMatrix(&backwardMatrix, N);
    SolveBackwardMatrix(L, U, backwardMatrix, P, N);
    printf("\n4) Other:\nA^-1:\n");
    printMatrix(backwardMatrix, N);

    //A * A^-1
    double** AA = nullptr;
    createMatrix(&AA, N);
    fillMatrixAsEmpty(AA, N);
    matrixMul(A, backwardMatrix, AA, N);
    printf("\nA * A^-1:\n");
    printMatrix(AA, N, true);
    double** AAminusE = nullptr;
    createMatrix(&AAminusE, N);
    matrixSub(AA, E, AAminusE, N);
    printf("(A * A^-1) - E:\n");
    printMatrix(AAminusE, N, true);

    //найти определитель
    auto det = computeDet(L, N, sign);
    printf("\ndet(A) = %f\n", det);

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

    printf("cubCond(A) = %f\noctCond(A) = %f\neuclidCond(A) = %f\n", cubCond, octCond, euclidCond);
    return 0;
}