#include <iostream>
#include <cmath>
#include <vector>
#include <cstring>
#include <omp.h>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <mpi.h>

using namespace std;

const double mineps = pow(10, -15);
const double RhoMin = pow(10, -60);
static bool info;

double globtime = 0.0;

struct Matrix
{
    int rows = 0;
    int cols = 0;
    int M = 0; // sizeA/sizeJA
    int N = 0; // sizeIA
    vector<int> IA = vector<int>();
    vector<int> JA = vector<int>();
    vector<double> A = vector<double >();

    Matrix(int rows, int cols, double *mtr = nullptr) : rows(rows), cols(cols), M(0)
    {

        N = max(cols, rows);
        if (mtr) {
            int m = 0;
            for (int i = 0; i < rows; i++) {
                bool first = true;
                for (int j = 0; j < cols; j++) {
                    if (mtr[i * cols + j]) {
                        if (first) {
                            first = false;
                            IA.push_back(m);
                        }
                        JA.push_back(j);
                        A.push_back(mtr[i * cols + j]);
                        m++;
                    }

                }
            }
            IA.push_back(m);
            M = m;
        } else {
            IA = vector<int>(rows + 1);
            JA = vector<int>();
            A = vector<double>();
        }
    }

    double operator()(int i, int j)
    {
        for (int k = IA[i]; k < IA[i + 1]; k++)
            if (j == JA[k])
                return A[k];
        return 0.0;
    }

    double get(int i, int j) const
    {
        for (int k = IA[i]; k < IA[i + 1]; k++)
            if (j == JA[k])
                return A[k];
        return 0.0;
    }


    void printCSR()
    {
        cout << "A" << endl;
        for (int i = 0; i < M; i++)
            cout << A[i] << ' ';
        cout << endl;

        cout << "JA" << endl;
        for (int i = 0; i < M; i++)
            cout << JA[i] << ' ';
        cout << endl;

        cout << "IA" << endl;
        for (int i = 0; i < N + 1; i++)
            cout << IA[i] << ' ';
        cout << endl;

        return;
    }

    void print(int flg = 0)
    {
        if (flg == 1)
            printCSR();

        cout << "Matrix" << endl;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                cout << get(i, j) << ' ';
            }
            cout << endl;
        }

        return;
    }

    Matrix& generator(int nx, int ny, int nz) {
        N = nx * ny * nz;
        int corner = 4 * 8; // 3 soseda
        int edge = 5 * 4 * ((nx - 2) + (ny - 2) + (nz - 2)); // 4 soseda
        int face = 6 * 2 * (((nx - 2) * (ny - 2)) + ((nx - 2) * (nz - 2)) + ((ny - 2) * (nz - 2))); // Xy Xz Yz
        int inner = 7 * (nx - 2) * (ny - 2) * (nz - 2); // 6 sosedey
        M = corner + edge + face + inner;

        this->IA = vector<int>(N + 1, 0);
        this->JA = vector<int>(M, 0);
        this->A = vector<double>(M, 0);
        IA[0] = 0;

        int ia = 0, ija = 0, iia = 1;
        int cnt = 0;

        for (int k = 0; k < nz; k++)
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++)
                {
                    int num = i + nx * j + nx * ny * k;
                    double sum = 0.0;
                    // double s = sin(i + j + 1);

                    if (k > 0)
                    {
                        // A[ia] = s;
                        A[ia] = sin(2 * num - nx * ny + 1);
                        sum += fabs(A[ia++]);
                        JA[ija++] = num - nx * ny;
                        cnt++;
                    }

                    if (j > 0)
                    {
                        // A[ia] = s;
                        A[ia] = sin(2 * num - nx + 1);
                        sum += fabs(A[ia++]);
                        JA[ija++] = num - nx;
                        cnt++;
                    }

                    if (i > 0)
                    {
                        // A[ia] = s;
                        A[ia] = sin(2 * num - 1 + 1);
                        sum += fabs(A[ia++]);
                        JA[ija++] = num - 1;
                        cnt++;
                    }

                    int tia = ia++, tija = ija++;

                    if (i < nx - 1)
                    {
                        // A[ia] = s;
                        A[ia] = sin(2 * num + 1 + 1);
                        sum += fabs(A[ia++]);
                        JA[ija++] = num + 1;
                        cnt++;
                    }

                    if (j < ny - 1)
                    {
                        // A[ia] = s;
                        A[ia] = sin(2 * num + nx + 1);
                        sum += fabs(A[ia++]);
                        JA[ija++] = num + nx;
                        cnt++;
                    }

                    if (k < nz - 1)
                    {
                        // A[ia] = s;
                        A[ia] = sin(2 * num + nx * ny + 1);
                        sum += fabs(A[ia++]);
                        JA[ija++] = num + nx * ny;
                        cnt++;
                    }

                    A[tia] = fabs(sum) * 1.1;
                    JA[tija] = num;
                    IA[iia++] = ++cnt;
                }
        return *this;
    }

    Matrix& diagonal(Matrix mat, int N) {

        for (int i = 0; i < mat.N; i++)
            if (mat.get(i, i) > 0.000001)
                M++;

        N = mat.N;

        IA = vector<int>(N + 1);
        JA = vector<int>(M);
        A = vector<double>(M);
        IA[0] = 0;
        int cnt = 0;

        int ia = 0, iia = 1;

        for (int i = 0; i < N; i++)
        {
            double diag = mat.get(i, i);

            if (diag > 0.000001)
            {
                A[ia] = 1.0 / diag;
                // A[ia] = diag;
                JA[ia++] = i;
                cnt++;
            }

            IA[iia++] = cnt;
        }
        return *this;
    }

    ~Matrix()
    {
        A.clear();
        JA.clear();
        IA.clear();

        N = 0;
        M = 0;
    }
};

double dot(const vector<double> &vec1, const vector<double> &vec2) {
    if (vec1.size() != vec2.size())
        return -1;
    double res = 0;

#pragma omp parallel for reduction(+:res)
    for (int i = 0; i < vec1.size(); i++) {
        res += vec1[i] * vec2[i];
    }
    return res;


}

int SpMV(Matrix &matrixA, const vector<double> &vectorB, vector<double> &vectorY) {
    int rowsA = matrixA.rows, colsA = matrixA.cols;
    auto size = vectorB.size();
    vector<double> vectorC(rowsA, 0);

    if (colsA != size)
    {
        cout << "Incorrect sizes!" << endl;
        return -1;
    }

#pragma omp parallel for
    for (int i = 0; i < matrixA.N; i++)
    {
        vectorC[i] = 0.0;
        double res = 0.0;
#pragma omp parallel for reduction(+:res)
        for (int j = matrixA.IA[i]; j < matrixA.IA[i + 1]; j++)
            res += matrixA.A[j] * vectorB[matrixA.JA[j]];
        vectorC[i] = res;
    }
    vectorY = vectorC;
    return 0;
}

void axpby(vector<double> &vectorX, const vector<double> &vectorY, double a, double b) {
    if (vectorX.size() != vectorY.size())
        return;
#pragma omp parallel for
    for (auto i = 0; i < vectorY.size(); i++){
        vectorX[i] = a * vectorX[i] + b * vectorY[i];
    }

}

struct result
{
    vector<double> res;
    int nit;

    result(vector<double> r, int i) {
        res = r;
        nit = i;
    }

    result() {
        res = {};
        nit = 0;
    }
};

int solver(const int N, Matrix &A, vector<double> &BB, int maxit, double tol, result &data) {
    vector<double> XX(N, 0);
    vector<double> RR = BB;
    vector<double> RR2 = BB;
    vector<double> PP, PP2, TT, VV, SS, SS2;
    double initres = sqrt(dot(BB, BB));
    double eps = max(mineps, tol*initres);
    double res = initres;
    double Rhoi_1 = 1.0, alphai = 1.0, wi = 1.0, betai_1 = 1.0, Rhoi_2 = 1.0, alphai_1 = 1.0, wi_1 = 1.0;
    Matrix DD(N, N);
    DD.diagonal(A, N);
    int I = 0;
//    printf("SOLVER_BICGSTAB: initres: %e; eps: %e; N = %d\n", initres, eps, N);
    for (I = 0; I < maxit; I++) {
        if (info) printf("It %d: res = %e tol=%e\n", I, res, res / initres);

        if (res < eps) {
            data = result(XX, I);
            break;
        }
        if (res > initres / mineps) return -1;

        if (I == 0) Rhoi_1 = initres*initres;
        else Rhoi_1 = dot(RR2, RR);
        if (fabs(Rhoi_1) < RhoMin) return -1;

        if (I == 0) PP = RR;
        else {
            betai_1 = (Rhoi_1 * alphai_1) / (Rhoi_2 * wi_1);
            axpby(PP, RR, betai_1, 1.0); // p=r+betai_1*(p-w1*v)
            axpby(PP, VV, 1.0, -1 * wi_1 * betai_1);
        }

        SpMV(DD, PP, PP2);

        SpMV(A, PP2, VV);

        alphai = dot(RR2, VV);
        if (fabs(alphai) < RhoMin) return -3;
        alphai = Rhoi_1 / alphai;

        SS = RR; // s=r-alphai*v
        axpby(SS, VV, 1.0, -1 * alphai);

        SpMV(DD, SS, SS2);

        SpMV(A, SS2, TT);

        wi = dot(TT, TT);
        if(fabs(wi) < RhoMin) return -4;
        wi = dot(TT, SS) / wi;
        if(fabs(wi) < RhoMin) return -5;

        // x=x+alphai*p2+wi*s2
        axpby(XX, PP2, 1.0, alphai);
        axpby(XX, SS2, 1.0, wi);

        RR = SS; // r=s-wi*t

        axpby(RR, TT, 1.0, -1 * wi);

        alphai_1 = alphai;
        Rhoi_2 = Rhoi_1;
        wi_1 = wi;

        res = sqrt(dot(RR, RR));
    }
    if (info) printf("Solver_BiCGSTAB: outres: %g tol: %g\n",res, tol);
    if (info) printf("Solver finished in %d iterations, res = %g tol = %e\n", I, res, tol);
    return I;
}

int main(int argc, char* argv[])
{
    if (argc < 8)
    {
        cout << "> Wrong number of in params, check your command" << endl;
        cout << "<Nx> <Ny> <Nz> <tol> <maxit> <omp> <Px> <Py> <Pz> <debug>" << endl;
        MPI_Finalize();
        return -1;
    }

    MPI_Init(&argc, &argv);

    int nProc = 0, myRank = 0;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);


    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int nz = atoi(argv[3]);
    double tol = atof(argv[4]);
    int maxit = atoi(argv[5]);
    int nt = atoi(argv[6]);
    bool qa = atoi(argv[7]);
    info = qa;
    long long int N = nx * ny * nz;

    Matrix test(N, N);
    test.generator(nx, ny, nz);
    vector<double> X(N), Y(N);
    for (int i = 0; i < N; i++) {
        X[i] = sin(i);
        Y[i] = cos(i);
    }
    vector<double> B(N);
    vector<double> Res(N);
    for (int i = 0; i < N; i++) {
        B[i] = sin(i);
    }
    result res;

    int Ntest = 10;
    vector<tuple<int, int, int>> mount = {make_tuple(1, 10, 100), make_tuple(1, 100, 100), make_tuple(10, 100, 100), make_tuple(100, 100, 100)};
    printf("testing sequential ops (N = ):\n");
    omp_set_num_threads(1);
    for (auto it : mount) {

       int nx = get<0>(it), ny = get<1>(it), nz = get<2>(it);
       int N = nx * ny * nz;

       double taxpyseq=0.0, tspmvseq=0.0, tdotseq=0.0, tsolverseq=0.0;
       const double axpyflop = Ntest*Ntest*N*3*1E-9;
       const double dotflop = Ntest*Ntest*N*2*1E-9;

       Matrix test(N, N);
       test.generator(nx, ny, nz);
       int M = test.A.size();
       const double spmvflop = (Ntest * Ntest)*M*2*1E-9;
       vector<double> X(N), Y(N);
       for (int i = 0; i < N; i++) {
           X[i] = sin(i);
           Y[i] = cos(i);
       }
       vector<double> B(N);
       vector<double> Res(N);
       for (int i = 0; i < N; i++) {
           B[i] = sin(i);
       }
       result res;
       omp_set_num_threads(1);
       for(int i=0; i<Ntest; i++){

           double t = omp_get_wtime();
           for(int j=0; j<Ntest; j++) dot(X, Y);
//            dot(X, Y);
           tdotseq += omp_get_wtime() - t;

           t = omp_get_wtime();
           for(int j=0; j<Ntest; j++) axpby(X, Y, 1.00001, 0.99999);
//            axpby(X, Y, 1.00001, 0.99999);
           taxpyseq += omp_get_wtime() - t;

           t = omp_get_wtime();
           for(int j=0; j<Ntest; j++) SpMV(test, B, Res);
//            SpMV(test, B, Res);
           tspmvseq += omp_get_wtime() - t;

       }
       printf("Sequential ops timing: (N = %d)\n", N);
       printf("dot time=%6.5fs GFLOPS=%6.2f\n", tdotseq / (Ntest * Ntest), dotflop/tdotseq);
       printf("axpy time=%6.5fs GFLOPS=%6.2f\n", taxpyseq / (Ntest * Ntest), axpyflop/taxpyseq);
       printf("SpMV time=%6.5fs GFLOPS=%6.2f\n", tspmvseq / (Ntest * Ntest), spmvflop/tspmvseq);

       double t = omp_get_wtime();
       solver(N, test, B, maxit, tol, res);
       tsolverseq = omp_get_wtime() - t;
       printf("solver time = %6.3fs\n", tsolverseq);

//        const int NTR = omp_get_num_procs();
       const int NTR = 10;
       for(int ntr=2; ntr<=NTR; ntr+=2) {
           for (int i = 0; i < N; i++) {
               X[i] = sin(i);
               Y[i] = cos(i);
           }
           for (int i = 0; i < N; i++) {
               B[i] = sin(i);
           }
           printf("testing parallel ops for ntr=%d:\n", ntr);
           omp_set_num_threads(ntr);
           double taxpypar=0.0, tspmvpar=0.0, tdotpar=0.0, tsolverpar=0.0;
           for(int i=0; i<Ntest; i++) {

               double t = omp_get_wtime();
               for(int j=0; j<Ntest; j++) dot(X, Y);
//                dot(X, Y);
               tdotpar += omp_get_wtime() - t;

               t = omp_get_wtime();
               for(int j=0; j<Ntest; j++) axpby(X, Y, 1.00001, 0.99999);
//                axpby(X, Y, 1.00001, 0.99999);
               taxpypar += omp_get_wtime() - t;

               t = omp_get_wtime();
               for(int j=0; j<Ntest; j++) SpMV(test, B, Res);
//                SpMV(test, B, Res);
               tspmvpar += omp_get_wtime() - t;

           }

           printf("dot time=%6.5fs GFLOPS=%6.2f Speedup=%6.2fX \n",
                  tdotpar / (Ntest * Ntest), dotflop/tdotpar, tdotseq/tdotpar);
           printf("axpy time=%6.5fs GFLOPS=%6.2f Speedup=%6.2fX \n",
                  taxpypar / (Ntest * Ntest), axpyflop/taxpypar, taxpyseq/taxpypar);
           printf("SpMV time=%6.5fs GFLOPS=%6.2f Speedup=%6.2fX \n",
                  tspmvpar / (Ntest * Ntest), spmvflop/tspmvpar, tspmvseq/tspmvpar);

           double t = omp_get_wtime();
           int r = solver(N, test, B, maxit, tol, res);
           tsolverpar = omp_get_wtime() - t;
           printf("solver time = %6.3fs\n", tsolverpar);
       }
       cout << endl;
   }

    MPI_Finalize();
    return 0;
}
