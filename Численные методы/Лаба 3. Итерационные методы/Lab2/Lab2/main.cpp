#include <iostream>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <iomanip>
#define M_PI            3.14159265358979323846

using namespace std;
const double minValue = 0.00001;
const double eps = 0.0001;

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

//копируем значения одной матрицы в другую
void copyVectorToVector(double* srcVector, double* dstVector, int N) {
    for (int i = 0; i < N; i++) {
        dstVector[i] = srcVector[i];
    }
}

//очистка вектора
void vectorClear(double* A, int N)
{
    for (int i = 0; i < N; i++)
        A[i] = 0.0;
}

//умножение матрицы на вектор
void matrixMulVec(double** A, double* B, double* C, int N)
{
    vectorClear(C, N);
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

//умножение матрицы на скаляр
void matrixMulScalar(double** A, double scalar, double** C, int N)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            C[i][j] = A[i][j] * scalar;
}

//вычитание матрицы
void matrixSub(double** A, double** B, double** C, int N, double scalar = 1.0)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            C[i][j] = A[i][j] - B[i][j] * scalar;
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

//сложение векторов
void vectorAdd(double* A, double* B, double* C, int N)
{
    for (int i = 0; i < N; i++)
        C[i] = A[i] + B[i];
}

//вычитание векторов
void vectorSub(double* A, double* B, double* C, int N)
{
    for (int i = 0; i < N; i++)
        C[i] = A[i] - B[i];
}

//умножение вектора на скаляр
void vectorMulScalar(double* A, double scalar, double* C, int N)
{
    for (int i = 0; i < N; i++)
        C[i] = A[i] * scalar;
}

//скалярное произведенеие векторов
double vectorMul(double* v1, double* v2, int N)
{
    double sum = 0;
    for (int i = 0; i < N; i++)
        sum += v1[i] * v2[i];
    return sum;
}

//норма вектора
double vectorEuNorm(double* vector, int N)
{
    double sum = 0;
    for (int i = 0; i < N; i++)
        sum += pow(vector[i], 2);
    return sqrt(sum);
}

double vectorNorm(double* v, int N)
{
    double max = abs(v[0]);
    for (int i = 1; i < N; i++)
        if (abs(v[i]) > max)
            max = abs(v[i]);
    return max;
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

//транспонировать матрицу
void transpose(double** matrix, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            swap(matrix[j][i], matrix[i][j]);
        }
    }
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

//чтение входных данных
void readInputData(const char* filename, double*** pA, double** pB, int* pN) {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::exception();
    int idx = 0;
    double value;
    file >> *pN;

    int N = *pN;
    createMatrix(pA, N);
    *pB = new double[N];
    while (file >> value)
    {
        if (idx < N * N) {
            (*pA)[idx / N][idx % N] = value;
        }
        else {
            (*pB)[idx % N] = value;
        }
        idx++;
    }
}

void printIterHeader()
{
    cout << setw(8) << " " << setw(8) << " " << setw(7) << " " << "|" << setw(12) << "Норма" << "|" << setw(12) << "Оценка" << "|\n";
    cout << setw(7) << "Itr" << "|" << setw(7) << "Tau" << "|" << setw(7) << "q" << "|" << setw(12) << "невязки" << "|" << setw(12) << "погрешности" << "|";
    cout << setw(12) << "x[1]" << "|" << setw(12) << "x[2]" << "|" << setw(12) << "x[3]" << "|" << setw(12) << "x[4]" << "|\n";
}

void printIterStep(int itr, double tau, double q, double norm, double err, double* v, int N, bool isAlpha = false, double alpha = 0.0)
{
    cout << setw(7) << itr << "|";
    cout << setw(7) << setprecision(4) << fixed << tau << "|";
    cout << setw(7) << setprecision(3) << fixed << q << "|";
    cout << setw(12) << setprecision(7) << fixed << norm << "|";
    cout << setw(12) << setprecision(7) << fixed << err << "|";
    for (int i = 0; i < N; i++)
        cout << setw(12) << setprecision(5) << fixed << v[i] << "|";
    if (isAlpha)
        cout << setw(12) << setprecision(7) << fixed << alpha;
    cout << endl;
}

//вычислить норму вектора-невязки
double calcVectorRNorm(double** A, double* B, double* X, double* prevX, double* tempVector, int N) {
    //R = Ax
    matrixMulVec(A, prevX, tempVector, N);
    //R = b - Ax
    vectorSub(B, tempVector, tempVector, N);
    //|R|
    return vectorNorm(tempVector, N) / vectorNorm(X, N);
}

//оценка нормы матрицы перехода q
double calcQ(double* X, double* prevX, double* tempVector, double* prevDelta, int N) {
    vectorSub(X, prevX, tempVector, N);
    auto delta = vectorNorm(tempVector, N);
    auto q = delta / *prevDelta;
    *prevDelta = delta;
    return q;
}

//считаем оценку погрешности
double calcError(double* X, double* prevX, double* tempVector, double q, int N)
{
    vectorSub(X, prevX, tempVector, N);
    return vectorNorm(tempVector, N) * (1 - q) / q;
}

//метод простой итерации
void SimpleIterationMethod(double** A, double normA, double* B, double* X, int N)
{
    auto prevX = new double[N];
    auto tempVector = new double[N];
    int itr = 1;
    double tau = 0.9 * 2 / normA;
    double normR = 0;
    auto prevDelta = 1.0;
    vectorClear(prevX, N);
    
    printIterHeader();
    do {
        //R = Ax
        matrixMulVec(A, prevX, tempVector, N);
        //R = b - Ax
        vectorSub(B, tempVector, tempVector, N);
        //(b - Ax) * tau
        vectorMulScalar(tempVector, tau, tempVector, N);
        //x(k+1) = x(k) + (b - Ax) * tau
        vectorAdd(prevX, tempVector, X, N);


        //норма невязки
        normR = calcVectorRNorm(A, B, X, prevX, tempVector, N);

        //оценка нормы матрицы перехода q
        auto q = calcQ(X, prevX, tempVector, &prevDelta, N);

        //погрешность
        auto err = calcError(X, prevX, tempVector, q, N);

        //вывод
        printIterStep(itr++, tau, q, normR, err, X, N);
        copyVectorToVector(X, prevX, N);
    } while (normR > eps);
}

//метод наискорейшего спуска
void FastDescentMethod(double** A, double* B, double* X, int N)
{
    auto prevX = new double[N];
    auto vectorR = new double[N];
    auto tempVector = new double[N];
    int itr = 1;
    double normR = 0;
    auto prevDelta = 1.0;
    vectorClear(prevX, N);

    printIterHeader();
    do {
        //R = Ax
        matrixMulVec(A, prevX, vectorR, N);
        //R = b - Ax
        vectorSub(B, vectorR, vectorR, N);
        //Ar
        matrixMulVec(A, vectorR, tempVector, N);
        //tau = <r, r> / <Ar, r>
        double tau = vectorMul(vectorR, vectorR, N) / vectorMul(tempVector, vectorR, N);
        //(b - Ax) * tau
        vectorMulScalar(vectorR, tau, tempVector, N);
        //x(k+1) = x(k) + (b - Ax) * tau
        vectorAdd(prevX, tempVector, X, N);


        //норма невязки
        normR = calcVectorRNorm(A, B, X, prevX, tempVector, N);

        //оценка нормы матрицы перехода q
        auto q = calcQ(X, prevX, tempVector, &prevDelta, N);

        //погрешность
        auto err = calcError(X, prevX, tempVector, q, N);

        //вывод
        printIterStep(itr++, tau, q, normR, err, X, N);
        copyVectorToVector(X, prevX, N);
    } while (normR > eps);
}

//считаем alpha
double calcAlpha(double tau, double prevTau, double prevAlpha, double R2, double prevR2)
{
    return 1.0 / (1.0 - R2 * tau / (prevR2 * prevTau * prevAlpha));
}

//метод сопряженных градиентов
void ConjugateGradientMethod(double** A, double* B, double* X, int N)
{
    auto prevPrevX = new double[N];
    auto prevX = new double[N];
    auto vectorR = new double[N];
    auto tempVector = new double[N];
    int itr = 1;
    double normR = 0;
    auto prevAlpha = 1.0;
    auto prevDelta = 1.0;
    auto prevR2 = 0.0;
    auto prevTau = 0.0;
    vectorClear(prevX, N);
    for (int i = 0; i < N; i++) prevPrevX[i] = 1.0;

    printIterHeader();
    do {
        //R = Ax
        matrixMulVec(A, prevX, vectorR, N);
        //R = b - Ax
        vectorSub(B, vectorR, vectorR, N);
        //Ar
        matrixMulVec(A, vectorR, tempVector, N);
        //tau = <r, r> / <Ar, r>
        double R2 = vectorMul(vectorR, vectorR, N);
        double tau = R2 / vectorMul(tempVector, vectorR, N);

        auto alpha = (itr == 1) ? 1.0 : calcAlpha(tau, prevTau, prevAlpha, R2, prevR2);
        prevR2 = R2;
        prevAlpha = alpha;
        prevTau = tau;

        //alpha * x(k)
        vectorMulScalar(prevX, alpha, X, N);
        //(1 - alpha) * x(k - 1)
        vectorMulScalar(prevPrevX, 1 - alpha, tempVector, N);
        vectorAdd(X, tempVector, X, N);
        //alpha * tau * R
        vectorMulScalar(vectorR, alpha * tau, tempVector, N);
        vectorAdd(X, tempVector, X, N);


        //норма невязки
        normR = calcVectorRNorm(A, B, X, prevX, tempVector, N);

        //оценка нормы матрицы перехода q
        auto q = calcQ(X, prevX, tempVector, &prevDelta, N);

        //погрешность
        auto err = calcError(X, prevX, tempVector, q, N);

        //вывод
        printIterStep(itr++, tau, q, normR, err, X, N, true, alpha);
        copyVectorToVector(prevX, prevPrevX, N);
        copyVectorToVector(X, prevX, N);
    } while (normR > eps);
}

int SOR(double** A, double* B, double* X, double omega, double N, double eps, bool printSteps);
double findBestOmega(double** A, double* B, int N) {
    auto X = new double[N];
    double omega = 0.1;
    double bestOmega = 0.1;
    int minItr = INT_MAX;
    while (omega < 2.0) {
        auto itr = SOR(A, B, X, omega, N, 0.01, false);
        if (itr < minItr)
        {
            minItr = itr;
            bestOmega = omega;
        }
        printf("w = %3.1f\tItr = %d\n", omega, itr);
        omega += 0.1;
    }

    printf("w* = %3.1f\tItrMin = %d\n", bestOmega, minItr);
    return bestOmega;
}

void SOR_calcVectorX(double** matrixA, double* vectorXNext, double* vectorX, double* vectorB, int N, double omega)
{
    for (int i = 0; i < N; i++)
    {
        double sumOne = 0;
        double sumTwo = 0;

        for (int j = 0; j < i; j++)
        {
            sumOne += matrixA[i][j] * vectorXNext[j];
        }
        for (int j = i + 1; j < N; j++)
        {
            sumTwo += matrixA[i][j] * vectorX[j];
        }

        double xWave = (vectorB[i] - sumOne - sumTwo) / matrixA[i][i];
        vectorXNext[i] = vectorX[i] + omega * (xWave - vectorX[i]);
    }
}

//метод ПВР
int SOR(double** A, double* B, double* X, double omega, double N, double eps, bool printSteps = true)
{
    auto prevX = new double[N];
    auto tempVector = new double[N];
    int itr = 1;
    double normR = 0;
    auto prevDelta = 1.0;
    vectorClear(prevX, N);

    if (printSteps) {
        printIterHeader();
    }

    do {
        //расчитываем текущий вектор X
        SOR_calcVectorX(A, X, prevX, B, N, omega);

        //норма невязки
        normR = calcVectorRNorm(A, B, X, prevX, tempVector, N);

        if (printSteps) {
            //оценка нормы матрицы перехода q
            auto q = calcQ(X, prevX, tempVector, &prevDelta, N);

            //погрешность
            auto err = calcError(X, prevX, tempVector, q, N);

            //вывод
            printIterStep(itr, omega, q, normR, err, X, N);
        }
        itr++;
        copyVectorToVector(X, prevX, N);
    } while (normR > eps);

    return itr;
}

void initSomeStuff(double** A, double* vectorB, double* vectorX, double** backwardMatrix, int N)
{
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
    int* P = new int[N]; //"матрица" перестановок (на самом деле подстановка)
    double sign = 1.0;

    LUdecomposition(L, U, P, rank, sign, N);

    if (rank != N)
        throw std::exception();

    SolveSOLE(L, U, vectorX, P, vectorB, N);
    SolveBackwardMatrix(L, U, backwardMatrix, P, N);
}

int main()
{
    system("chcp 1251");

    //чтение
    double** A = nullptr;
    double* B = nullptr;
    int N;
    readInputData("in2.txt", &A, &B, &N);
    auto X_LU = new double[N];
    auto X_Method1 = new double[N];
    auto X_Method2 = new double[N];
    auto X_Method3 = new double[N];
    auto X_Method4 = new double[N];

    //обратная матрица
    double** backwardMatrix = nullptr;
    createMatrix(&backwardMatrix, N);

    //LU, ...
    initSomeStuff(A, B, X_LU, backwardMatrix, N);

    //транспонируем матрицу
    double** trA = nullptr;
    createMatrix(&trA, N);
    copyMatrixToMatrix(A, trA, N);
    transpose(trA, N);
    double** trBackA = nullptr;
    createMatrix(&trBackA, N);
    copyMatrixToMatrix(backwardMatrix, trBackA, N);
    transpose(trBackA, N);

    printf("Variant=\nb\n");
    printVector(B, N);
    printf("A\n");
    printMatrix(A, N);
    double euclidNorm = computEuclidNorm(A, trA, N);
    double euclidNorm2 = computEuclidNorm(backwardMatrix, trBackA, N);
    printf("Норма матрицы=%f\n", euclidNorm);

    printf("Метод простой итерации\n");
    SimpleIterationMethod(A, euclidNorm, B, X_Method1, N);

    printf("\nМетод наискорейшего спуска\n");
    FastDescentMethod(A, B, X_Method2, N);

    printf("\nМетод ПВР - выбор оптимального w\n");
    auto bestOmega = findBestOmega(A, B, N);

    printf("\nМетод ПВР\n");
    SOR(A, B, X_Method3, bestOmega, N, eps);

    printf("\nМетод сопряженных градиентов\n");
    ConjugateGradientMethod(A, B, X_Method4, N);

    auto condA = euclidNorm * euclidNorm2;
    cout << "\nЧисло обусловленности " << condA << endl;
    cout << "\nТеоретическая оценка чиcла итераций ";
    cout << "\nМетод простых итераций " << (long)(log(1 / eps) / 2 * condA);
    cout << "\nМетод наискорейшего спуска " << (long)(log(1 / eps) * condA);
    cout << "\nМетод ПВР " << (long)(log(1 / eps) / 4 * sqrt(condA));
    cout << "\nМетод сопряженных градиентов " << (long)(log(2 / eps) / 2 * sqrt(condA));

    auto tempVector = new double[N];
    cout << "\n\nСравнение с LU разложением\n";
    printVector(X_LU, 0);
    cout << "\nРазница между LU разложением и методом простых итераций\n";
    vectorSub(X_LU, X_Method1, tempVector, N);
    printVector(tempVector, N, true);
    cout << "\nРазница между LU разложением и методом наискорейшего спуска\n";
    vectorSub(X_LU, X_Method2, tempVector, N);
    printVector(tempVector, N, true);
    cout << "\nРазница между LU разложением и методом ПВР\n";
    vectorSub(X_LU, X_Method3, tempVector, N);
    printVector(tempVector, N, true);
    cout << "\nРазница между LU разложением и методом сопряженных градиентов\n";
    vectorSub(X_LU, X_Method4, tempVector, N);
    printVector(tempVector, N, true);
    return 0;
}