#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <time.h> 
using namespace std;

const int N = 4;
const double gamma = 0.9;
const double epsilon = 0.0001;
const char* fileInput = "testin.txt"; // "testin.txt"; //
const char* fileOutput = "output2.txt";//"testout.txt"; //
ofstream output;
double condA = 0;
double normA = 0;

double** readMatrix() //считывание исходной матрицы
{
    ifstream file;
    file.open(fileInput, ios::in);
    double** m = new double* [N];
    for (int i = 0; i < N; i++)
        m[i] = new double[N];
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            file >> m[i][j];
    file.close();
    return m;
}
void printMatrix(double** m, bool exp) // печать матрицы в файл
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            if (exp) output << setw(20) << scientific << m[i][j];
            else output << setw(12) << setprecision(7) << fixed << m[i][j];
        output << endl;
    }

}

void printVector(double* v, bool exp) // печать вектора в файл
{
    for (int j = 0; j < N; j++)
        if (exp) output << setw(20) << scientific << v[j];
        else output << setw(12) << setprecision(7) << fixed << v[j];
    output << endl;
}
void printHead()
{
    output << setw(8) << " " << setw(8) << " " << setw(7) << " " << "|" << setw(12) << "Норма" << "|" << setw(12) << "Оценка" << "|\n";
    output << setw(7) << "Itr" << "|" << setw(7) << "Tau" << "|" << setw(7) << "q" << "|" << setw(12) << "невязки" << "|" << setw(12) << "погрешности" << "|";
    output << setw(12) << "x[1]" << "|" << setw(12) << "x[2]" << "|" << setw(12) << "x[3]" << "|" << setw(12) << "x[4]" << "|\n";
}

void printRow(int itr, double tau, double q, double norm, double err, double* v, bool isAlpha, double alpha)
{
    output << setw(7) << itr << "|";
    output << setw(7) << setprecision(4) << fixed << tau << "|";
    output << setw(7) << setprecision(3) << fixed << q << "|";
    output << setw(12) << setprecision(7) << fixed << norm << "|";
    output << setw(12) << setprecision(7) << fixed << err << "|";
    for (int i = 0; i < N; i++)
        output << setw(12) << setprecision(5) << fixed << v[i] << "|";
    if (isAlpha)
        output << setw(12) << setprecision(7) << fixed << alpha;
    output << endl;
}
double** createAndCopyMatrix(double** A) //создание аналогичной матрицы
{
    double** B = new double* [N];
    for (int i = 0; i < N; i++)
        B[i] = new double[N];
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            B[i][j] = A[i][j];
    return B;
}
void deleteMatrix(double** m) //удаление матрицы, очищение памяти
{
    for (int i = 0; i < N; i++)
        delete[] m[i];
    delete[] m;
}

double** multMatrix(double** A, double** B) // перемножение двух матриц
{
    double** res = new double* [N];
    for (int i = 0; i < N; i++)
        res[i] = new double[N];

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            res[i][j] = 0;
            for (int k = 0; k < N; k++)
                res[i][j] += A[i][k] * B[k][j];
        }
    return res;
}
double* multMatrVec(double** A, double* b) //  перемножения матрицы и вектора
{
    double* v = new double[N];
    for (int i = 0; i < N; i++)
    {
        v[i] = 0;
        for (int k = 0; k < N; k++)
            v[i] += A[i][k] * b[k];
    }
    return v;
}

double** subMatrices(double** a, double** b)
{
    double** res = new double* [N];
    for (int i = 0; i < N; i++)
        res[i] = new double[N];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            res[i][j] = a[i][j] - b[i][j];
        }
    }
    return res;
}
double* subVec(double* a, double* b)
{
    double* v = new double[N];
    for (int i = 0; i < N; i++)
        v[i] = a[i] - b[i];
    return v;
}
// скалярное произведение векторов
double scalarMult(double* a, double* b)
{
    double sum = 0;
    for (int i = 0; i < N; i++)
        sum += a[i] * b[i];
    return sum;
}

void transposeMatrix(double** m) // транспонированная матрица
{
    for (int i = 0; i < N - 1; i++)
        for (int j = i + 1; j < N; j++)
            swap(m[i][j], m[j][i]);

}
double EuclideanNorm(double** A)  //Евклидова норма
{
    double max = 0, elem_max = 0, fi;
    int k = 0, i_max, j_max;
    bool f = true;
    double** At = createAndCopyMatrix(A);
    transposeMatrix(At);
    double** AtA = multMatrix(At, A);
    double** H = new double* [N];
    double** Ht = new double* [N];;
    for (int i = 0; i < N; i++)
    {
        H[i] = new double[N];
        Ht[i] = new double[N];
    }

    while (f)
    {
        elem_max = 0;
        for (int i = 0; i < N - 1; i++)
            for (int j = i + 1; j < N; j++)
                if (abs(AtA[i][j]) > elem_max)
                {
                    elem_max = abs(AtA[i][j]);
                    i_max = i;
                    j_max = j;
                }
        f = abs(elem_max) > 0.0001;
        if (f)
        {
            double temp = 2 * AtA[i_max][j_max] / (AtA[i_max][i_max] - AtA[j_max][j_max]);
            fi = 0.5 * atan(temp);
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    H[i][j] = i == j ? 1 : 0;
            H[i_max][i_max] = cos(fi);
            H[j_max][j_max] = cos(fi);
            H[i_max][j_max] = -sin(fi);
            H[j_max][i_max] = sin(fi);
            deleteMatrix(Ht);
            Ht = createAndCopyMatrix(H);
            transposeMatrix(Ht);
            double** AtA2 = multMatrix(Ht, AtA);
            double** AtA3 = multMatrix(AtA2, H);
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    AtA[i][j] = AtA3[i][j];
            deleteMatrix(AtA2);
            deleteMatrix(AtA3);


        }
    }

    for (int i = 0; i < N; i++)
        if (abs(AtA[i][i]) > max)
            max = abs(AtA[i][i]);

    deleteMatrix(At);
    deleteMatrix(AtA);
    deleteMatrix(H);
    deleteMatrix(Ht);
    return sqrt(max);
}
double EuclidianVector(double* b)
{
    double res = 0;
    for (int i = 0; i < N; i++) {
        res += b[i] * b[i];
    }
    return sqrt(res);
}
double normVector(double* v)
{
    double max = abs(v[0]);
    for (int i = 1; i < N; i++)
        if (abs(v[i]) > max)
            max = abs(v[i]);
    return max;
}
double* solveSystem(double** U, double** L, double* b)
{
    double y[N];
    for (int i = 0; i < N; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= y[j] * L[i][j];
        }
        y[i] /= L[i][i];
    }

    double* x_new = new double[N];
    for (int i = N - 1; i >= 0; i--) {
        x_new[i] = y[i];
        for (int j = i + 1; j < N; j++) {
            x_new[i] -= U[i][j] * x_new[j];
        }

    }
    return x_new;
}
double** inverseMatrix(double** U, double** L, double** P) //обратная матрица
{
    double** A_inverse = new double* [N];
    for (int i = 0; i < N; i++)
        A_inverse[i] = new double[N];
    double* e = new double[N];
    double* Pe;
    for (int k = 0; k < N; k++)
    {

        for (int i = 0; i < N; i++)
            e[i] = 0;
        e[k] = 1;
        Pe = multMatrVec(P, e);
        double* x = solveSystem(U, L, Pe);

        for (int i = 0; i < N; i++)
            A_inverse[i][k] = x[i];
        delete[]x;
        delete[]Pe;
    }
    delete[]e;
    return A_inverse;

}

void swapRows(double** m, int a, int b) // перестановка строк
{
    double temp;
    for (int i = 0; i < N; i++)
    {
        temp = m[a][i];
        m[a][i] = m[b][i];
        m[b][i] = temp;
    }
}
int indexMax(double** m, int j) //номер строки с макс. элементом по столбцу
{
    int inx = j;
    for (int i = j; i < N; i++)
        if (abs(m[i][j]) > abs(m[inx][j]))
            inx = i;
    return inx;
}

void getLU(double** U, double** L, double** P, double** A)
{
    double** subA = createAndCopyMatrix(A);
    double maxE;
    int i, j, k, m;
    int r = N; // ранг матрицы
    for (k = 0; k < N; k++)
    {

        m = indexMax(U, k); //номер строки с макс. элементом
        if (k != m) //поменять строки местами
        {
            swapRows(U, m, k);
            swapRows(L, m, k);
            swapRows(P, m, k);
            swapRows(subA, m, k);
            //countswaps++;
            output << " k = " << k + 1 << " m = " << m + 1 << endl;
        }


        maxE = U[k][k];
        for (j = k; j < N; j++)
            U[k][j] /= maxE;
        if (abs(maxE) < 0.00001) r = k;
        for (i = k + 1; i < N; i++)
        {
            double tmp = U[i][k];
            for (int j = k; j < N; j++)
                U[i][j] = U[i][j] - U[k][j] * tmp;

        }
        for (int j = 0; j <= k; j++) {
            L[k][j] = subA[k][j];
            for (int i = 0; i < j; i++)
                L[k][j] -= L[k][i] * U[i][j];
        }
        /*output << " U = " << endl;
        printMatrix(U, 0);
        output << " L = " << endl;
        printMatrix(L, 0);*/
    }
    //output << " Rang(A) = " << r << endl;
    deleteMatrix(subA);
}
//значение нормы невязки
double getNormResidual(double** A, double* x_next, double* x, double* b)
{
    return normVector(subVec(multMatrVec(A, x), b)) / normVector(x_next);

}
//считаем оценку нормы матрицы перехода q;
double getNormQ(double* x_next, double* x, double* x_prev)
{
    return normVector(subVec(x_next, x)) / normVector(subVec(x, x_prev));
}
// считаем оценку погрешности
double getErr(double* x_next, double* x, double q)
{
    return normVector(subVec(x_next, x)) * (1 - q) / q;

}
void copyVectors(double* x_next, double* x, double* x_prev)
{
    for (int i = 0; i < N; i++)
    {
        x_prev[i] = x[i];
        x[i] = x_next[i];
    }
}

double getTau(double** A, double* x, double* b)
{
    double* r = subVec(multMatrVec(A, x), b);
    double* Ar = multMatrVec(A, r);
    return scalarMult(r, r) / scalarMult(Ar, r);
}

double* methodSimpleIterations(double** A, double* b)
{
    output << "\n Метод простой итерации \n";
    printHead();
    double* x = new double[N] { 0.0, 0.0, 0.0, 0.0 };
    double x_prev[N] = { 1.0, 1.0, 1.0 ,1.0 };
    double x_next[N] = { 0.0, 0.0, 0.0 ,0.0 };
    int itr = 1;
    double q, norm, err;
    double tau = gamma * 2 / normA;
    do
    {
        //получение следующего значения вектора
        for (int i = 0; i < N; i++)
        {
            double temp = 0;
            x_next[i] = 0;
            for (int j = 0; j < N; j++)
                temp += A[i][j] * x[j];
            x_next[i] = x[i] + tau * (b[i] - temp);
        }
        //норма невязки
        norm = getNormResidual(A, x_next, x, b);

        //норма вектора q
        q = getNormQ(x_next, x, x_prev);

        //погрешность
        err = getErr(x_next, x, q);
        printRow(itr, tau, q, norm, err, x_next, 0, 0);
        itr++;
        copyVectors(x_next, x, x_prev);

    } while (norm > epsilon);
    return x;
}

double* methodSteepestDescent(double** A, double* b)
{
    output << "\n Метод наискорейшего спуска\n";
    printHead();
    double* x = new double[N] { 0.0, 0.0, 0.0, 0.0 };
    double x_prev[N] = { 1.0, 1.0, 1.0 ,1.0 };
    double x_next[N] = { 0.0, 0.0, 0.0 ,0.0 };
    double* r;
    int itr = 1;
    double q, norm, err;
    double tau;
    do
    {
        tau = getTau(A, x, b);

        r = subVec(multMatrVec(A, x), b);

        //получение следующего значения вектора
        for (int i = 0; i < N; i++)
        {
            double temp = 0;
            x_next[i] = 0;
            for (int j = 0; j < N; j++)
                temp += A[i][j] * x[j];
            x_next[i] = x[i] + tau * (b[i] - temp);
            /* x_next[i] = x[i] - tau * r[i];*/
        }

        //норма невязки
        norm = getNormResidual(A, x_next, x, b);

        //норма вектора q
        q = getNormQ(x_next, x, x_prev);

        //погрешность
        err = getErr(x_next, x, q);

        printRow(itr, tau, q, norm, err, x_next, 0, 0);
        itr++;
        copyVectors(x_next, x, x_prev);
        delete[] r;
    } while (norm > epsilon);


    return x;
}

void getNextVectSOR(double** A, double* x_next, double* x, double* b, double omega)
{
    //формула с лекции 
    for (int i = 0; i < N; i++)
    {
        double sumOne = 0;
        double sumTwo = 0;

        for (int j = 0; j < i; j++)
            sumOne += A[i][j] * x_next[j];
        for (int k = i + 1; k < N; k++)
            sumTwo += A[i][k] * x[k];

        double x_next_wave = (1 / A[i][i]) * (b[i] - sumOne - sumTwo);

        x_next[i] = x[i] + omega * (x_next_wave - x[i]);
    }
}
double getOptimalOmega(double** A, double* b)
{
    output << "\n Метод ПВР - выбор оптимального w\n";
    double omega = 0.1;
    double optimalPEPEGA = omega;
    double eps = 0.01;
    double norm = 0;
    int minItr = INT_MAX;

    while (omega < 2)
    {
        double x[N] = { 0.0, 0.0, 0.0, 0.0 }; //k
        double x_next[N] = { 0.0, 0.0, 0.0 ,0.0 }; //k+1
        int itr = 0;

        do {
            itr++;
            getNextVectSOR(A, x_next, x, b, omega);
            norm = getNormResidual(A, x_next, x, b);

            for (int i = 0; i < N; i++)
                x[i] = x_next[i];

        } while (norm > eps);

        if (itr < minItr)
        {
            minItr = itr;
            optimalPEPEGA = omega;
        }

        output << " w=" << omega << " Itr=" << itr << endl;
        omega += 0.1;
    }

    output << endl;
    output << " w*= " << optimalPEPEGA << " ItrMin=" << minItr << endl;
    return optimalPEPEGA;
}
double* methodSOR(double** A, double* b)
{
    double omega = getOptimalOmega(A, b);
    double* x = new double[N] { 0.0, 0.0, 0.0, 0.0 }; //k
    double x_prev[N] = { 1.0, 1.0, 1.0 ,1.0 }; //k-1
    double x_next[N] = { 0.0, 0.0, 0.0 ,0.0 }; //k+1
    int itr = 0;
    double norm, q, err;

    output << "\n Метод ПВР\n";
    printHead();

    do
    {
        itr++;

        //получение следующего значения вектора
        getNextVectSOR(A, x_next, x, b, omega);

        //норма вектора q
        q = getNormQ(x_next, x, x_prev);

        //погрешность
        err = getErr(x_next, x, q);

        //норма невязки
        norm = getNormResidual(A, x_next, x, b);

        copyVectors(x_next, x, x_prev);

        printRow(itr, omega, q, norm, err, x_next, 0, 0);

    } while (norm > epsilon);

    return x;
}

double getScalarR(double** A, double* x, double* b)
{
    double* r = subVec(multMatrVec(A, x), b);
    return scalarMult(r, r);
}
double* methodConjugateGradient(double** A, double* b)
{
    double* x = new double[N] { 0.0, 0.0, 0.0, 0.0 }; //k
    double x_prev[N] = { 1.0, 1.0, 1.0 ,1.0 }; //k-1
    double x_next[N] = { 0.0, 0.0, 0.0 ,0.0 }; //k+1
    int itr = 0;
    double norm, q, err, tau, tauPrev, alpha = 1;
    output << "\n Метод сопряженных градиентов\n";
    printHead();

    double scalarR;
    double scalarRPrev;
    double* r;

    do {
        itr++;
        tau = getTau(A, x, b);
        scalarR = getScalarR(A, x, b);
        scalarRPrev = getScalarR(A, x_prev, b);
        r = subVec(multMatrVec(A, x), b);

        //на первом шаге alpha = 1. Формулы по лекции(в соотв. топике)
        if (itr != 1)
            alpha = 1 / (1 - (tau * scalarR) / (tauPrev * alpha * scalarRPrev));

        for (int i = 0; i < N; i++)
            x_next[i] = alpha * x[i] + (1 - alpha) * x_prev[i] - tau * alpha * r[i];

        //норма вектора q
        q = getNormQ(x_next, x, x_prev);
        //оценка погрешности
        err = getErr(x_next, x, q);
        //норма невязки
        norm = getNormResidual(A, x_next, x, b);
        copyVectors(x_next, x, x_prev);
        tauPrev = tau;
        printRow(itr, tau, q, norm, err, x_next, 1, alpha);

    } while (norm > epsilon);

    return x;
}
int main()
{
    srand(time(0));
    output.open(fileOutput);
    double** A = readMatrix();              //исходная матрица А
    double** U = createAndCopyMatrix(A);    // матрица U изначально равна исходной
    double** L = new double* [N];           // матрица L
    double** P = new double* [N];           //матрица перестановок
    double** A_reverse;                     // А^(-1)
    double b[] = { 1, 2, 3, 4 };

    for (int i = 0; i < N; i++)
    {
        P[i] = new double[N];
        L[i] = new double[N];
    }
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            L[i][j] = 0;                    // матрица L изначально нулевая
            P[i][j] = i == j ? 1 : 0;       // матрица P изначально единичная
        }
    getLU(U, L, P, A);
    double* Pb = multMatrVec(P, b);
    double* x = solveSystem(U, L, Pb);
    A_reverse = inverseMatrix(U, L, P);

    output << " Variant = 27b" << endl;
    output << " b:" << endl;
    printVector(b, 0);
    output << " A:" << endl;
    printMatrix(A, 0);
    normA = EuclideanNorm(A);
    condA = normA * EuclideanNorm(A_reverse);
    output << "\n Норма матрицы = " << normA << endl;

    double* x_simpleIter = methodSimpleIterations(A, b);
    double* x_descent = methodSteepestDescent(A, b);
    double* x_SOR = methodSOR(A, b);
    double* x_Gradient = methodConjugateGradient(A, b);


    output << "\n Число обусловленности = " << setw(10) << setprecision(5) << fixed << condA << endl;
    output << "\n Теоретическая оценка чиcла итераций ";
    output << "\n Метод простых итераций " << int(log(1 / epsilon) / 2 * condA);
    output << "\n Метод наискорейшего спуска " << int(log(1 / epsilon) * condA);
    output << "\n Метод ПВР " << int(log(1 / epsilon) / 4 * sqrt(condA));
    output << "\n Метод сопряженных градиентов " << int(log(2 / epsilon) / 2 * sqrt(condA));

    output << "\n\n LU разложение \n";
    printVector(x, 0);
    output << "\n Разница между LU разложением и методом простых итераций \n";
    printVector(subVec(x, x_simpleIter), 0);
    output << "\n Разница между LU разложением и методом наискорейшего спуска \n";
    printVector(subVec(x, x_descent), 0);
    output << "\n Разница между LU разложением и методом ПВР \n";
    printVector(subVec(x, x_SOR), 0);
    output << "\n Разница между LU разложением и методом сопряженных градиентов \n";
    printVector(subVec(x, x_Gradient), 0);


    deleteMatrix(A);
    deleteMatrix(U);
    deleteMatrix(L);
    deleteMatrix(P);
    delete[] Pb;
    delete[] x;
    delete[] x_simpleIter;
    delete[] x_descent;
    delete[] x_SOR;
    delete[] x_Gradient;
    output.close();
}

