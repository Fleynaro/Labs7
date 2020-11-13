#include <iostream>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <functional>
#include <iomanip>
#include <cstdarg>

using namespace std;

//Вектор
class Vector
{
    double* m_values;
public:
    int m_N;

    Vector(double* values, int N)
        : m_values(values), m_N(N)
    {}

    double& operator[](int index) {
        return m_values[index];
    }

    void fillBy(double value = 0) {
        for (int i = 0; i < m_N; i++) {
            m_values[i] = value;
        }
    }
};

//Матрица
class Matrix
{
    double* m_values;
public:
    int m_N;
    int m_M;

    Matrix(double* values, int N, int M)
        : m_values(values), m_N(N), m_M(M)
    {}

    Matrix(int N)
        : m_N(N), m_M(N)
    {
        m_values = new double[N * N];
        clear();
    }

    Vector operator[](int index) {
        return Vector(&m_values[m_N * index], m_N);
    }

    void print() {
        printf("\n");
        for (int i = 0; i < m_N; i++) {
            for (int j = 0; j < m_M; j++) {
                printf("%.7f ", (*this)[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    void clear() {
        for (int i = 0; i < m_M; i++) {
            (*this)[i].fillBy(0);
        }
    }

    static void Add(Matrix& C, Matrix& A, Matrix& B, double sign = 1.0) {
        for (int i = 0; i < C.m_N; i++) {
            for (int j = 0; j < C.m_M; j++) {
                C[i][j] = A[i][j] + B[i][j] * sign;
            }
        }
    }

    static void Mul(Matrix& C, Matrix& A, Matrix& B) {
        for (int i = 0; i < C.m_N; i++) {
            for (int j = 0; j < C.m_M; j++) {
                C[i][j] = 0.0;
                for (int k = 0; k < A.m_M; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
    }

    template<int N>
    static Matrix Create(double values[N][N]) {
        Matrix matrix(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                matrix[i][j] = values[i][j];
            }
        }
        return matrix;
    }
};



/*
    Алгоритм Штрассена
    1) На входе: 2 квадратных матрицы A и B (в виде одномерного массива) На выходе: матрица C
    2) Надо использовать минимум памяти и её аллокации
    3) Обозначить границу n, начиная с которой идет переключение на обычный алгоритм
    4) Выделить объем памяти для каждой итерации для 3-х матриц a, b, c. Пусть изначальный размер массива n
        - первая итерация потребует 3 * (n / 2) ^ 2
        - вторая итерация потребует 3 * (n / 4) ^ 2
        - ... до тех пор, пока не будет достигнута отсечка: n / (2^h) <= m

    На каждой итерации:
    1. обнуляем матрицу C
    2. используем 3 матрицы размера 1/4 текущей: A, B, C
    3. к обнуленной матрице C прибавляем только
    4. входные 2 матрицы не должны меняться
*/

class StrassenAlgorithm
{
    Matrix* m_A;
    Matrix* m_B;
    Matrix* m_C;
    int m_switchBoundary;
    double* m_tempData;
public:
    StrassenAlgorithm(Matrix* A, Matrix* B, Matrix* C, int switchBoundary = 2)
        : m_A(A), m_B(B), m_C(C), m_switchBoundary(switchBoundary)
    {
        allocateMemoryForTempData();
    }

    void start() {
        recMul(*m_A, *m_B, *m_C, 0);
    }
private:
    void recMul(Matrix& A, Matrix& B, Matrix& C, int depth) {
        auto matrixDstOperand = getTempMatrixByIndex(0, depth);
        if(matrixDstOperand.m_N < m_switchBoundary) {
            Matrix::Mul(C, A, B);
            return;
        }

        C.clear();
        auto matrixSrcOperand1 = getTempMatrixByIndex(1, depth);
        auto matrixSrcOperand2 = getTempMatrixByIndex(2, depth);

        //P1
        takeQuarter(matrixSrcOperand1, A, 0, 0);
        takeQuarter(matrixSrcOperand2, A, 1, 1);
        
        Matrix::Add(matrixSrcOperand1, matrixSrcOperand1, matrixSrcOperand2, 1.0);
        takeQuarter(matrixDstOperand, B, 0, 0);
        takeQuarter(matrixSrcOperand2, B, 1, 1);
        Matrix::Add(matrixSrcOperand2, matrixDstOperand, matrixSrcOperand2, 1.0);

        recMul(matrixSrcOperand1, matrixSrcOperand2, matrixDstOperand, depth + 1);
        addQuarter(C, matrixDstOperand, 0, 0);
        addQuarter(C, matrixDstOperand, 1, 1);
        

        //P2
        takeQuarter(matrixSrcOperand1, A, 1, 0);
        takeQuarter(matrixSrcOperand2, A, 1, 1);
        Matrix::Add(matrixSrcOperand1, matrixSrcOperand1, matrixSrcOperand2, 1.0);
        takeQuarter(matrixSrcOperand2, B, 0, 0);
        
        recMul(matrixSrcOperand1, matrixSrcOperand2, matrixDstOperand, depth + 1);
        addQuarter(C, matrixDstOperand, 1, 0);
        addQuarter(C, matrixDstOperand, 1, 1, -1.0);


        //P3
        takeQuarter(matrixSrcOperand1, B, 0, 1);
        takeQuarter(matrixSrcOperand2, B, 1, 1);
        Matrix::Add(matrixSrcOperand2, matrixSrcOperand1, matrixSrcOperand2, -1.0);
        takeQuarter(matrixSrcOperand1, A, 0, 0);

        recMul(matrixSrcOperand1, matrixSrcOperand2, matrixDstOperand, depth + 1);
        addQuarter(C, matrixDstOperand, 0, 1);
        addQuarter(C, matrixDstOperand, 1, 1);


        //P4
        takeQuarter(matrixSrcOperand1, B, 1, 0);
        takeQuarter(matrixSrcOperand2, B, 0, 0);
        Matrix::Add(matrixSrcOperand2, matrixSrcOperand1, matrixSrcOperand2, -1.0);
        takeQuarter(matrixSrcOperand1, A, 1, 1);

        recMul(matrixSrcOperand1, matrixSrcOperand2, matrixDstOperand, depth + 1);
        addQuarter(C, matrixDstOperand, 0, 0);
        addQuarter(C, matrixDstOperand, 1, 0);


        //P5
        takeQuarter(matrixSrcOperand1, A, 0, 0);
        takeQuarter(matrixSrcOperand2, A, 0, 1);
        Matrix::Add(matrixSrcOperand1, matrixSrcOperand1, matrixSrcOperand2, 1.0);
        takeQuarter(matrixSrcOperand2, B, 1, 1);

        recMul(matrixSrcOperand1, matrixSrcOperand2, matrixDstOperand, depth + 1);
        addQuarter(C, matrixDstOperand, 0, 0, -1.0);
        addQuarter(C, matrixDstOperand, 0, 1);


        //P6
        takeQuarter(matrixSrcOperand1, A, 1, 0);
        takeQuarter(matrixSrcOperand2, A, 0, 0);
        Matrix::Add(matrixSrcOperand1, matrixSrcOperand1, matrixSrcOperand2, -1.0);
        takeQuarter(matrixDstOperand, B, 0, 0);
        takeQuarter(matrixSrcOperand2, B, 0, 1);
        Matrix::Add(matrixSrcOperand2, matrixDstOperand, matrixSrcOperand2, 1.0);

        recMul(matrixSrcOperand1, matrixSrcOperand2, matrixDstOperand, depth + 1);
        addQuarter(C, matrixDstOperand, 1, 1);


        //P7
        takeQuarter(matrixSrcOperand1, A, 0, 1);
        takeQuarter(matrixSrcOperand2, A, 1, 1);
        Matrix::Add(matrixSrcOperand1, matrixSrcOperand1, matrixSrcOperand2, -1.0);
        takeQuarter(matrixDstOperand, B, 1, 0);
        takeQuarter(matrixSrcOperand2, B, 1, 1);
        Matrix::Add(matrixSrcOperand2, matrixDstOperand, matrixSrcOperand2, 1.0);

        recMul(matrixSrcOperand1, matrixSrcOperand2, matrixDstOperand, depth + 1);
        addQuarter(C, matrixDstOperand, 0, 0);
    }

    void takeQuarter(Matrix& dst, Matrix& src, int row, int col) {
        for (int i = 0; i < dst.m_N; i++) {
            for (int j = 0; j < dst.m_N; j++) {
                dst[i][j] = src[row * dst.m_N + i][col * dst.m_N + j];
            }
        }
    }

    void addQuarter(Matrix& dst, Matrix& src, int row, int col, double sign = 1.0) {
        for (int i = 0; i < src.m_N; i++) {
            for (int j = 0; j < src.m_N; j++) {
                dst[row * src.m_N + i][col * src.m_N + j] += src[i][j] * sign;
            }
        }
    }

    Matrix getTempMatrixByIndex(int matrixIdx, int depth) {
        int offset = 0;
        int n = getSourceMatrixSize();
        int k = 2;
        while (depth > 0) {
            offset += 3 * (n * n) / (k * k);
            k *= 2;
            depth--;
        }
        auto matrixSize = n / k;
        return Matrix(&m_tempData[offset + matrixIdx * matrixSize * matrixSize], matrixSize, matrixSize);
    }

    void allocateMemoryForTempData() {
        int size = 0;
        int n = getSourceMatrixSize();
        int k = 2;
        while (n / k >= m_switchBoundary) {
            size += 3 * (n * n) / (k * k);
            k *= 2;
        }
        m_tempData = new double[size];
    }

    int getSourceMatrixSize() {
        return m_A->m_N;
    }
};

int main()
{
    system("chcp 1251");

    double valuesA[8][8] = {
        { 1.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 2.0 },
        { 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0 },
        { 20.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 2.0 },
        { 3.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 2.0 },
        { 3.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 2.0 },
        { 3.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 2.0 },
        { 3.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 2.0 },
        { 3.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 2.0 }
    };
    auto A = Matrix::Create<8>(valuesA);
    //for (int i = 0; i < A.m_N; i++) for (int j = 0; j < A.m_N; j++) A[i][j] = (i == j) ? 1.0 : 0.0;
    A.print();

    double valuesB[8][8] = {
        { 1.0, 3.0, 0.0, 2.0, 1.0, 0.0, 0.0, 2.0 },
        { 0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0 },
        { 20.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 2.0 },
        { 3.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 2.0 },
        { 3.0, 0.0, 0.0, 1.0, 1.0, 3.0, 0.0, 3.0 },
        { 3.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 2.0 },
        { 3.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 2.0 },
        { 3.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 3.0 }
    };
    auto B = Matrix::Create<8>(valuesB);
    //for (int i = 0; i < B.m_N; i++) for (int j = 0; j < B.m_N; j++) B[i][j] = (i == j) ? 1.0 : 0.0;
    B.print();

    Matrix C1(8);
    Matrix C2(8);

    Matrix::Mul(C1, A, B);

    //todo: привидение матрицы к нужной степени
    StrassenAlgorithm strassenAlgorithm(&A, &B, &C2, 2);
    strassenAlgorithm.start();
    
    printf("\nResult:");
    C1.print();
    C2.print();
}