#include <iostream>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <sstream>
#include <functional>
#include <iomanip>
#include <cstdarg>
#include <mutex>
#include <atomic>
#include <list>
#include <Windows.h>

using namespace std;

enum class Algorithm {
    Simple = 1,
    Strassen
};

//Вектор
class Vector
{
    double* m_values;
    double m_value = 0.0;
public:
    int m_N;

    Vector(double* values, int N)
        : m_values(values), m_N(N)
    {}

    virtual double& operator[](int index) {
        return (index < m_N) ? m_values[index] : m_value;
    }

    void fillBy(double value = 0) {
        for (int i = 0; i < m_N; i++) {
            m_values[i] = value;
        }
    }
};

class AbstractMatrix
{
public:
    virtual int getN() = 0;

    virtual int getM() = 0;

    virtual Vector operator[](int index) = 0;

    void clear() {
        for (int i = 0; i < getM(); i++) {
            (*this)[i].fillBy(0);
        }
    }

    void print() {
        printf("\n");
        for (int i = 0; i < getN(); i++) {
            for (int j = 0; j < getM(); j++) {
                printf("%.7f ", (*this)[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    static void Add(AbstractMatrix& C, AbstractMatrix& A, AbstractMatrix& B, double sign = 1.0) {
        for (int i = 0; i < C.getN(); i++) {
            for (int j = 0; j < C.getM(); j++) {
                C[i][j] = A[i][j] + B[i][j] * sign;
            }
        }
    }

    static void Mul(AbstractMatrix& C, AbstractMatrix& A, AbstractMatrix& B) {
        for (int i = 0; i < C.getN(); i++) {
            for (int j = 0; j < C.getM(); j++) {
                C[i][j] = 0.0;
                for (int k = 0; k < A.getM(); k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
    }

    static void Equal(AbstractMatrix& A, AbstractMatrix& B) {
        for (int i = 0; i < A.getN(); i++) {
            for (int j = 0; j < A.getM(); j++) {
                if (A[i][j] != B[i][j])
                    throw std::logic_error("not equal!");
            }
        }
    }
};

//Матрица
class Matrix : public AbstractMatrix
{
    double* m_values;
    int m_N;
    int m_M;
public:
    Matrix(double* values, int N, int M)
        : m_values(values), m_N(N), m_M(M)
    {}

    Matrix(int N)
        : m_N(N), m_M(N)
    {
        m_values = new double[N * N];
        clear();
    }

    int getN() override {
        return m_N;
    }

    int getM() override {
        return m_M;
    }

    Vector operator[](int index) override {
        return Vector(&m_values[m_N * index], m_N);
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

    static void DestroyMatrix(Matrix* matrix) {
        delete[] matrix->m_values;
    }
};

//Расширенная матрица
class ExtMatrix : public AbstractMatrix
{
    AbstractMatrix* m_srcMatrix = nullptr;
    int m_extN = 0;
public:
    ExtMatrix() = default;

    ExtMatrix(AbstractMatrix* srcMatrix, int extN)
        : m_srcMatrix(srcMatrix), m_extN(extN)
    {}

    int getN() override {
        return m_extN;
    }

    int getM() override {
        return m_extN;
    }

    Vector operator[](int index) override {
        if (index >= m_srcMatrix->getN())
            return Vector(nullptr, 0);
        return (*m_srcMatrix)[index];
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
    ExtMatrix m_A;
    ExtMatrix m_B;
    ExtMatrix m_C;
    int m_switchBoundary;
    double* m_tempData;
public:
    StrassenAlgorithm(AbstractMatrix* A, AbstractMatrix* B, AbstractMatrix* C, int switchBoundary = 2)
        : m_switchBoundary(switchBoundary)
    {
        auto requiredN = GetNearestDegreeOf2(max(max(A->getN(), B->getN()), C->getN()));
        m_A = ExtMatrix(A, requiredN);
        m_B = ExtMatrix(B, requiredN);
        m_C = ExtMatrix(C, requiredN);
        allocateMemoryForTempData();
    }

    ~StrassenAlgorithm() {
        delete[] m_tempData;
    }

    void start() {
        recMul(m_A, m_B, m_C, 0);
    }
private:
    void recMul(AbstractMatrix& A, AbstractMatrix& B, AbstractMatrix& C, int depth) {
        auto matrixDstOperand = getTempMatrixByIndex(0, depth);
        if(matrixDstOperand.getN() < m_switchBoundary) {
            AbstractMatrix::Mul(C, A, B);
            return;
        }

        C.clear();
        auto matrixSrcOperand1 = getTempMatrixByIndex(1, depth);
        auto matrixSrcOperand2 = getTempMatrixByIndex(2, depth);

        //P1
        takeQuarter(matrixSrcOperand1, A, 0, 0);
        takeQuarter(matrixSrcOperand2, A, 1, 1);
        
        AbstractMatrix::Add(matrixSrcOperand1, matrixSrcOperand1, matrixSrcOperand2, 1.0);
        takeQuarter(matrixDstOperand, B, 0, 0);
        takeQuarter(matrixSrcOperand2, B, 1, 1);
        AbstractMatrix::Add(matrixSrcOperand2, matrixDstOperand, matrixSrcOperand2, 1.0);

        recMul(matrixSrcOperand1, matrixSrcOperand2, matrixDstOperand, depth + 1);
        addQuarter(C, matrixDstOperand, 0, 0);
        addQuarter(C, matrixDstOperand, 1, 1);
        

        //P2
        takeQuarter(matrixSrcOperand1, A, 1, 0);
        takeQuarter(matrixSrcOperand2, A, 1, 1);
        AbstractMatrix::Add(matrixSrcOperand1, matrixSrcOperand1, matrixSrcOperand2, 1.0);
        takeQuarter(matrixSrcOperand2, B, 0, 0);
        
        recMul(matrixSrcOperand1, matrixSrcOperand2, matrixDstOperand, depth + 1);
        addQuarter(C, matrixDstOperand, 1, 0);
        addQuarter(C, matrixDstOperand, 1, 1, -1.0);


        //P3
        takeQuarter(matrixSrcOperand1, B, 0, 1);
        takeQuarter(matrixSrcOperand2, B, 1, 1);
        AbstractMatrix::Add(matrixSrcOperand2, matrixSrcOperand1, matrixSrcOperand2, -1.0);
        takeQuarter(matrixSrcOperand1, A, 0, 0);

        recMul(matrixSrcOperand1, matrixSrcOperand2, matrixDstOperand, depth + 1);
        addQuarter(C, matrixDstOperand, 0, 1);
        addQuarter(C, matrixDstOperand, 1, 1);


        //P4
        takeQuarter(matrixSrcOperand1, B, 1, 0);
        takeQuarter(matrixSrcOperand2, B, 0, 0);
        AbstractMatrix::Add(matrixSrcOperand2, matrixSrcOperand1, matrixSrcOperand2, -1.0);
        takeQuarter(matrixSrcOperand1, A, 1, 1);

        recMul(matrixSrcOperand1, matrixSrcOperand2, matrixDstOperand, depth + 1);
        addQuarter(C, matrixDstOperand, 0, 0);
        addQuarter(C, matrixDstOperand, 1, 0);


        //P5
        takeQuarter(matrixSrcOperand1, A, 0, 0);
        takeQuarter(matrixSrcOperand2, A, 0, 1);
        AbstractMatrix::Add(matrixSrcOperand1, matrixSrcOperand1, matrixSrcOperand2, 1.0);
        takeQuarter(matrixSrcOperand2, B, 1, 1);

        recMul(matrixSrcOperand1, matrixSrcOperand2, matrixDstOperand, depth + 1);
        addQuarter(C, matrixDstOperand, 0, 0, -1.0);
        addQuarter(C, matrixDstOperand, 0, 1);


        //P6
        takeQuarter(matrixSrcOperand1, A, 1, 0);
        takeQuarter(matrixSrcOperand2, A, 0, 0);
        AbstractMatrix::Add(matrixSrcOperand1, matrixSrcOperand1, matrixSrcOperand2, -1.0);
        takeQuarter(matrixDstOperand, B, 0, 0);
        takeQuarter(matrixSrcOperand2, B, 0, 1);
        AbstractMatrix::Add(matrixSrcOperand2, matrixDstOperand, matrixSrcOperand2, 1.0);

        recMul(matrixSrcOperand1, matrixSrcOperand2, matrixDstOperand, depth + 1);
        addQuarter(C, matrixDstOperand, 1, 1);


        //P7
        takeQuarter(matrixSrcOperand1, A, 0, 1);
        takeQuarter(matrixSrcOperand2, A, 1, 1);
        AbstractMatrix::Add(matrixSrcOperand1, matrixSrcOperand1, matrixSrcOperand2, -1.0);
        takeQuarter(matrixDstOperand, B, 1, 0);
        takeQuarter(matrixSrcOperand2, B, 1, 1);
        AbstractMatrix::Add(matrixSrcOperand2, matrixDstOperand, matrixSrcOperand2, 1.0);

        recMul(matrixSrcOperand1, matrixSrcOperand2, matrixDstOperand, depth + 1);
        addQuarter(C, matrixDstOperand, 0, 0);
    }

    void takeQuarter(AbstractMatrix& dst, AbstractMatrix& src, int row, int col) {
        for (int i = 0; i < dst.getN(); i++) {
            for (int j = 0; j < dst.getN(); j++) {
                dst[i][j] = src[row * dst.getN() + i][col * dst.getN() + j];
            }
        }
    }

    void addQuarter(AbstractMatrix& dst, AbstractMatrix& src, int row, int col, double sign = 1.0) {
        for (int i = 0; i < src.getN(); i++) {
            for (int j = 0; j < src.getN(); j++) {
                dst[row * src.getN() + i][col * src.getN() + j] += src[i][j] * sign;
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
        return m_A.getN();
    }

    static int GetNearestDegreeOf2(int n) {
        n--;
        for (int i = 1; i <= 32; i++) {
            if ((n >>= 1) == 0)
                return 1 << i;
        }
        return 0;
    }
};

//Парсер
class MatrixParser
{
    std::ifstream m_file;
    int m_firstMatrixSize = -1;
public:
    MatrixParser(std::string filename = "in.txt")
        : m_file(filename)
    {}

    ~MatrixParser() {
        m_file.close();
    }

    Matrix parseMatrix() {
        if (!m_file.is_open())
            throw std::exception();
        int N = m_firstMatrixSize;
        if (N == -1) {
            m_file >> N;
            m_firstMatrixSize = N;
        }
        Matrix matrix(N);

        double value;
        int idx = 0;
        while (idx < N * N && m_file >> value)
        {
            matrix[idx / N][idx % N] = value;
            idx++;
        }

        return matrix;
    }
};

//Измерить время теста
struct Measure
{
    std::chrono::steady_clock::time_point m_start;

    void start() {
        m_start = std::chrono::high_resolution_clock::now();
    }

    long long elapsed() {
        auto elapsed = std::chrono::high_resolution_clock::now() - m_start;
        return std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    }
};

//Логгер тестов (формат .csv)
class TestLoggerCSV
{
    std::ofstream m_file;
    std::mutex m_mutex;
public:
    TestLoggerCSV(std::string filename = "tests.csv")
        : m_file(filename)
    {
        if (!m_file.is_open())
            throw std::exception();
        m_file.clear();
        m_file << "V" << "," << "testId" << "," << "algorithmId" << "," << "boundary" << "," << "time" << endl;
    }

    ~TestLoggerCSV() {
        m_file.close();
    }

    void insert(int V, int testId, Algorithm algorithmId, int strassenBoundary, Measure& measure) {
        auto elapsed = measure.elapsed();
        m_mutex.lock();
        std::stringstream ss;
        ss << V << "," << testId << "," << (int)algorithmId << "," << (int)strassenBoundary << "," << elapsed << endl;
        m_file << ss.str();
        cout << ss.str();
        m_mutex.unlock();
    }
};

class TestUnit
{
    TestLoggerCSV* m_testLoggerCSV;
    int m_testCount;
    int m_V;
    Algorithm m_algorithm;
    int m_strassenBoundary;
public:
    TestUnit(TestLoggerCSV* testLoggerCSV, int testCount, int V, Algorithm algorithmId, int strassenBoundary)
        : m_testLoggerCSV(testLoggerCSV), m_testCount(testCount), m_V(V), m_algorithm(algorithmId), m_strassenBoundary(strassenBoundary)
    {}

    void start() {
        for (int i = 0; i < m_testCount; i++) {
            std::thread task([&](int testId) {
                startOneTest(testId);
                }, i);
            task.join();
        }

    }

    static void GenerateMatrix(Matrix& A) {
        for (int i = 0; i < A.getN(); i++) {
            for (int j = 0; j < A.getM(); j++) {
                A[i][j] = rand() % 1000;
            }
        }
    }
private:

    void startOneTest(int testId) {
        Measure measure;
        Matrix A(m_V);
        Matrix B(m_V);
        Matrix C(m_V);
        GenerateMatrix(A);
        GenerateMatrix(B);

        if (m_algorithm == Algorithm::Simple) {
            measure.start();
            AbstractMatrix::Mul(C, A, B);
        }
        else {
            measure.start();
            StrassenAlgorithm strassenAlgorithm(&A, &B, &C, 64);
            strassenAlgorithm.start();
        }

        m_testLoggerCSV->insert(m_V, testId, m_algorithm, m_strassenBoundary, measure);
    }
};

class Test
{
    TestLoggerCSV* m_testLoggerCSV;
    std::atomic<int> m_completedTestCount;
    std::list<std::thread> m_tasks;
public:
    Test(TestLoggerCSV* testLoggerCSV)
        : m_testLoggerCSV(testLoggerCSV)
    {}

    void start()
    {
        //compareBothAlgorithms();
        compareBothAlgorithms2();

        //compareBothAlgorithmsOnSameVolumes();

        //checkOptimalParameterOfStrassen();

        for (auto& task : m_tasks)
            task.detach();

        wait();
    }

    //сравнить результат 2-х алгоритмов
    void compareBothAlgorithms2() {
        m_tasks.push_back(std::thread([&]() {
            task(10, 50, Algorithm::Simple);
            task(5, 250, Algorithm::Simple);
            task(5, 300, Algorithm::Simple);
            task(5, 400, Algorithm::Simple);
            task(5, 450, Algorithm::Simple);
            task(1, 550, Algorithm::Simple);
            task(1, 1024, Algorithm::Simple);
            m_completedTestCount++;
            }));
        m_tasks.push_back(std::thread([&]() {
            task(10, 50, Algorithm::Strassen, 16);
            task(5, 250, Algorithm::Strassen, 64);
            task(5, 300, Algorithm::Strassen, 64);
            task(5, 400, Algorithm::Strassen, 64);
            task(5, 450, Algorithm::Strassen, 64);
            task(1, 550, Algorithm::Strassen, 128);
            task(1, 1024, Algorithm::Strassen, 128);
            m_completedTestCount++;
            }));
    }

    //сравнить результат 2-х алгоритмов
    void compareBothAlgorithms() {
        m_tasks.push_back(std::thread([&]() {
            task(8, 50, Algorithm::Simple);
            task(8, 50, Algorithm::Strassen, 16);
            m_completedTestCount++;
            }));
        m_tasks.push_back(std::thread([&]() {
            task(5, 250, Algorithm::Simple);
            task(5, 250, Algorithm::Strassen, 64);
            m_completedTestCount++;
            }));
        m_tasks.push_back(std::thread([&]() {
            task(1, 512, Algorithm::Simple);
            task(1, 512, Algorithm::Strassen, 128);
            m_completedTestCount++;
            }));
        m_tasks.push_back(std::thread([&]() {
            task(1, 1024, Algorithm::Simple);
            task(1, 1024, Algorithm::Strassen, 128);
            m_completedTestCount++;
            }));
    }

    //Сравнить оба алгоритма на объеме 250
    void compareBothAlgorithmsOnSameVolumes() {
        m_tasks.push_back(std::thread([&]() {
            task(3, 460, Algorithm::Simple);
            //task(10, 250, Algorithm::Strassen, 64);
            m_completedTestCount++;
            }));

        m_tasks.push_back(std::thread([&]() {
            //task(10, 250, Algorithm::Simple);
            task(3, 460, Algorithm::Strassen, 64);
            m_completedTestCount++;
            }));
    }

    //проверить, какая отсечка лучше всего
    void checkOptimalParameterOfStrassen() {
        m_tasks.push_back(std::thread([&]() {
            task(5, 512, Algorithm::Strassen, 32);
            m_completedTestCount++;
            }));
        m_tasks.push_back(std::thread([&]() {
            task(5, 512, Algorithm::Strassen, 64);
            m_completedTestCount++;
            }));
        m_tasks.push_back(std::thread([&]() {
            task(5, 512, Algorithm::Strassen, 128);
            m_completedTestCount++;
            }));
    }

private:
    void task(int testCount, int V, Algorithm algorithmId, int strassenBoundary = 0) {
        TestUnit testUnit(m_testLoggerCSV, testCount, V, algorithmId, strassenBoundary);
        testUnit.start();     
    }

    bool isCompleted() {
        return m_completedTestCount == m_tasks.size();
    }

    void wait() {
        while (!isCompleted())
            Sleep(100);
    }
};

//равен ли результат алгоритмов
void checkCorrectness() {
    Matrix A(500);
    Matrix B(500);
    Matrix C1(500);
    Matrix C2(500);
    TestUnit::GenerateMatrix(A);
    TestUnit::GenerateMatrix(B);

    AbstractMatrix::Mul(C1, A, B);
    StrassenAlgorithm strassenAlgorithm(&A, &B, &C2, 64);
    strassenAlgorithm.start();

    //C1[144][244] = 0.2;
    AbstractMatrix::Equal(C1, C2);
}

int main()
{
    system("chcp 1251");

    TestLoggerCSV testLoggerCSV("matrix_mul_tests.csv");

    MatrixParser matrixParser("in.txt");
    Matrix A = matrixParser.parseMatrix();
    Matrix B = matrixParser.parseMatrix();
    Matrix C1(A.getN());
    Matrix C2(B.getN());

    printf("A:\n");
    A.print();
    printf("B:\n");
    B.print();

    Matrix::Mul(C1, A, B);
    StrassenAlgorithm strassenAlgorithm(&A, &B, &C2, 2);
    strassenAlgorithm.start();
    
    printf("\nResult:\n");
    printf("Тривиальный:\n");
    C1.print();
    printf("Оптимизированный:\n");
    C2.print();

    //checkCorrectness();


    Test test(&testLoggerCSV);
    test.start();
}