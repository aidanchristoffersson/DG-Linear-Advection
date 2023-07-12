#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

#include "DiscGalerkin.h"

// constant global variables
static const double a {0.5};
static const double xMin {0};
static const double xMax {40};

// sharpness of bump
static const double c {250};

// analytical solution // initial condition when t = 0
double func(const double x, const double t)
//{ return sin(x - (a * t)); }
{ return 0.5 * (1 + tanh(c * (x - (a * t) - 20))); }


// affine mapping
double affineMap(const double& r, const double& xL, const double& xR)
{
    return ((r + 1) / 2) * (xR - xL) + xL;
}

double mapAffine(const double& x, const double& xL, const double& xR)
{
    return 2 * ((x - xL) / (xR - xL)) - 1;
}

//template <size_t K> // size_t Np> // k x Np grid
class Grid {
private: // variables
    const unsigned int K;
    const unsigned int Nv;
    const unsigned int Np;
    const double h;

    std::unique_ptr<double[]> gridFlattened;
    std::unique_ptr<double *[]> grid; // each row is the set of nodes within each cell

public:
    std::unique_ptr<double[]> vX;

private: // functions
    void LogGrid() {
        for (unsigned int i = 0; i < this->K; i++) {
            for (unsigned int j = 0; j < Np; j++) {
                std::cout << grid[i][j] << " ";
            }
            std::cout << " - cell: " << i + 1 << std::endl;
        }
    }

public:
    // Constructor || input: # of cells, # of nodes per cell, boundaries
    Grid(const size_t &K, const size_t &Np, const double *glNodes)

            : K(K), Nv(K + 1), Np(Np), vX(new double[Nv]), h((xMax - xMin) / this->K),
              gridFlattened(std::make_unique<double[]>(this->K * this->Np)),
              grid(std::make_unique<double *[]>(this->K))

    {
        for (unsigned int i = 0; i < this->Nv; i++)
        {
            vX[i] = (h * i) + xMin; // assigning cell vertices (equally spaced)
        }

        for (unsigned int i = 0; i < this->K; i++) // iterating through rows
        {
            grid[i] = &(gridFlattened[i * Np]); // assigning indices s.t. a step is a jump to the next row

            for (unsigned int j = 0; j < Np; j++) {
                grid[i][j] = affineMap(glNodes[j], vX[i], vX[i + 1]);
            }
        }
//        LogGrid();
    }

    // indexing operator to access grid values
    double *&operator[](size_t index) { return grid[index]; }

    double *&operator[](size_t index) const { return grid[index]; }

    void write(const std::string &classification) const {
        // write grid
        std::ofstream x("Data/xGrid, " + classification + ".txt");

        if (!x.is_open()) { std::cout << "Error: Could not open file." << std::endl; }
        else {
            for (size_t i = 0; i < K * Np; i++) {
                x << gridFlattened[i] << " ";
            }
        }

        x.close();
    }

    void writeGridSpacing(const std::string &classification) const
    {
        // residual files
        std::ofstream dx("Data/GridSpacing, " + classification + ".txt", std::ios_base::app);

        if (!dx.is_open()) { std::cout << "Error: Could not open GridSpacing." << std::endl; }
        else {
            dx << h << " ";
        }

        dx.close();
    }

    [[nodiscard]] double dx() const
    { return h; }

};


class Solution
{
private: // variables
    const unsigned int K;
    const unsigned int Np;

    std::unique_ptr<double[]> uFlattened;
    std::unique_ptr<double*[]> u;

    std::unique_ptr<double[]> uAnalytical;

    double time {};
    double residual {};

    // may not be needed, only returning final time and residual value
    std::vector<double> residuals; //

private:

    [[maybe_unused]] void LogSoln()
    {
        for (unsigned int i = 0; i < this->K; i++) // same data structure as for the grid to have 1-1 mapping
        {
            for (unsigned int j = 0; j < Np; j++)
            {
                std::cout << u[i][j] << "-";
                std::cout << uAnalytical[(Np * i + j)] << " ";
            }
            std::cout << " - cell: " << i+1 << std::endl;
        }
    }

public:

    // constructor for initial condition
    Solution(const Grid& grid, const unsigned int& K, const unsigned int& N, const unsigned int& Np)

    : K(K), Np(Np), uFlattened(std::make_unique<double[]>(this->K * this->Np)),
    u(std::make_unique<double*[]>(this->K)), uAnalytical(std::make_unique<double[]>(this->K * this->Np))
    {
        for (unsigned int i = 0; i < this->K; i++) // same data structure as for the grid to have 1-1 mapping
        {
            u[i] = &(uFlattened[i * Np]);

            for (unsigned int j = 0; j < Np; j++) {
                // assign initial condition for both analytical and numerical solution

                u[i][j] = func(grid[i][j], (double)0);

                uAnalytical[(i * Np + j)] = u[i][j];
            }
        }
    }

    // indexing operator to access solution
    double*& operator[](size_t index) { return u[index]; }
    double*& operator[](size_t index) const { return u[index]; }

    double copy(size_t i, size_t j)
    { return u[i][j]; }

    double& Time()
    { return time; }

    void timestep(const double& dt)
    { time += dt; }

    // update analytical solution
    void analyticalSoln(const Grid& grid)
    {
        for (unsigned int i = 0; i < K; i++)
        {
            for (unsigned int j = 0; j < Np; j++)
            uAnalytical[(Np * i) + j] = func(grid[i][j], time);
        }
    }

    [[nodiscard]] double analytical(size_t i, size_t j) const  // input indices
    {
        return uAnalytical[(Np * i) + j];       // final index: Np * (K - 1) + (K - 1) = Np * K - Np + K - 1
    }

    void residualCalc() // append residual value
    {
        double resid {};

        for (unsigned int i = 0; i < K * Np; i++)
        {
            resid += pow(abs(uFlattened[i] - uAnalytical[i]), 2);  // rms error
        }
        residual = sqrt(resid / (K * Np));
        residuals.emplace_back(residual); // append rms error to residual vector
    }

    template <size_t Np>
    double TrueResidualCalc(const Grid& grid, const double* gLnodes, const LagrangianPoly<Np>& Lagrange) const
    {
        // number of points for plotting per cell (multiply by K for total points)
//        const size_t num = round((500 * 250) / K);
        const size_t num = 800;

        // initialize residual value
        double resid {};

        // placeholder variables
        double x;
        double r;

        // vertices
        double xL;
        double xR;

        // initialize step variable.
        double step;

        for (size_t i = 0; i < K; i++) // iterate over each cell
        {
            xL = grid.vX[i];
            xR = grid.vX[i + 1];

            step = (xR - xL) / (num - 1); // update step variable

            x = xL; // set x to leftmost point in cell

            for (size_t j = 0; j < num; j++) // sample at (# = num) points
            {
                r = mapAffine(x, xL, xR);

                // calculate residual at point
                resid += pow(Lagrange.L(r, gLnodes, u[i]) - func(x, time), 2);  // rms error

                // increment x each time
                x += step;
            }
        }
        return sqrt(resid / ((double)num * (double)K));
    }

    // writing solution to file
    void writeSoln(const std::string& classification, int count) const
    {
//        std::cout << "Time: " << time;
        // recording time

        // write time to file
        std::ofstream Time("Data/Times, " + classification + ".txt", std::ios_base::app);
        Time << time << " ";
        Time.close();

        // write iteration to file
        std::ofstream Count("Data/Count, " + classification + ".txt", std::ios_base::app);
        Count << count << " ";
        Count.close();

        std::string label = std::to_string(count);

        // solution files
        std::ofstream soln("Data/Solution, " + classification + ".txt");

        std::ofstream analyticalSoln("Data/AnalyticalSoln, " + classification + ".txt");

        if (!soln.is_open() || !analyticalSoln.is_open())
        { std::cout << "Error: Could not open file." << std::endl; }
        else
        {
            for (size_t i = 0; i < K * Np; i++)
            {
                soln << uFlattened[i] << " ";
                analyticalSoln << uAnalytical[i] << " ";
            }
        }
        soln.close();
        analyticalSoln.close();
    }

    template <size_t Np>
    void writeResidual(const std::string& classification, const Grid& grid, const double* gLnodes, const LagrangianPoly<Np>& Lagrange)
    {
//        // residual files
        std::ofstream convergence("Data/Residuals, " + classification + ".txt", std::ios_base::app);

        residual = TrueResidualCalc<Np>(grid, gLnodes, Lagrange);

        if (!convergence.is_open()) { std::cout << "Error: Could not open Residual." << std::endl; }
        else {
            convergence << residual << " ";
            }

        convergence.close();

        std::cout << "\nResidual at " << time << "s: "<< residual << std::endl;
    }
};

template <size_t K, size_t Np>
class AdvectionRHS {
private: // variables

    double numericalFLux[K + 1];

    double J;

    double rhsuFlattened[K * Np];
    double *rhsu[K];

private: // functions

    void NumericalFlux(Solution& u, const double* gLnodes, const LagrangianPoly<Np>& Lagrange)
    {
        // boundary domain
        double rR {1};
        double rL {-1} ;

        // boundary values
        double uL {};
        double uR {};

        // for leftmost vertex use internal flux of first element
//        numericalFLux[0] = - a * Lagrange.L(rL, gLnodes, u[0]);     // a * u(rL) * n(normal vector)

        numericalFLux[0] = a * u[0][0]; // Gauss Lobatto Points

        for (unsigned int i = 1; i < K; i++) // iterate through each internal boundary
        {
            // Gauss Legendre Calculation
//            uL = Lagrange.L(rR, gLnodes, u[i - 1]); // uL = u(rR)[cell: k]
//            uR = Lagrange.L(rL, gLnodes, u[i]);     // uR = u(rL)[cell: k + 1]

            // Gauss Lobatto Calculation
            uL = u[i - 1][Np - 1];
            uR = u[i][0];

            numericalFLux[i] = (a * (uL + uR) / 2) + (abs(a) * (uL - uR) / 2); // a > 0 -> a * uL, a < 0 -> a * uR
        }
        // update last flux (counteracting introduced error as using last node instead of boundary)
        // Gauss Legendre
//        uAend = u[K - 1][Np - 1];
//        numericalFLux[K] = a * uAend;

        numericalFLux[K] = a * u[K - 1][Np - 1];
//        std::cout << numericalFLux[K + 1] << std::endl;
    };


    void NumericalFluxK(double**& p, const double* gLnodes, const LagrangianPoly<Np>& Lagrange)
    {
        // due to Gauss-Lobatto, many inputs no longer required

        // boundary domain
        const double rR {1};
        const double rL {-1} ;

        // boundary values
        double uL {};
        double uR {};

        // for leftmost vertex use internal flux of first element
//        numericalFLux[0] = -a * Lagrange.L(rR, gLnodes, p[0]);

        numericalFLux[0] = a * p[0][0]; // Gauss Lobatto Points

        for (unsigned int i = 1; i < K; i++) // iterate through each internal boundary
        {
            // Gauss Legendre Calculation
//            uL = Lagrange.L(rR, gLnodes, p[i - 1]); // uL = u(rR)[cell: k]
//            uR = Lagrange.L(rL, gLnodes, p[i]);     // uR = u(rL)[cell: k + 1]

            // Gauss Lobatto Calculation
            uL = p[i - 1][Np - 1];
            uR = p[i][0];

            numericalFLux[i] = (a * (uL + uR) / 2) + (abs(a) * (uL - uR) / 2); // a > 0: a * uL, a < 0: a * uR
        }
//        // update last flux to counteract error introduced by fixing final element using GL nodes
//        numericalFLux[K] = a * uAend; // uAnalytical(xFinal, time) updated at beginning of each timestep

        // Gauss Lobatto end node
        numericalFLux[K] = a * p[K - 1][Np - 1];
    };

public:
    // constructor
    AdvectionRHS(Solution& u, const Grid& grid, const storedGLobatto<Np - 1>& gL, const LagrangianPoly<Np>& LagrangePoly)
    : J(grid.dx() / 2)
    {
        for (unsigned int i = 0; i < K; i++) // iterate over cells and assign row pointers
        {
            rhsu[i] = &(rhsuFlattened[i * Np]);
        }

        update(u, gL, LagrangePoly);
    }

    // update Numerical Flux and Advection RHS
    void update(Solution& u, const storedGLobatto<Np - 1> & gL, const LagrangianPoly<Np>& D)
    {
        NumericalFlux(u, gL.nodes, D);

        double Sf {};        // (S)^{T} * (au)
        double fLR {};       // f * lq @ xR (within cell k)
        double fLL {};       // f * lq @ xL (within cell k)


        for (unsigned int k = 0; k < K; k++) // iterate over cells
        {
            for (unsigned int q = 0; q < Np; q++) // iterate over lagrangian functions (row of vector)
            {
                Sf = 0; // placeholder value for each row

                // calculate stiffness * f vector (row)
                for (unsigned int p = 0; p < Np; p++) // iterate over nodes
                {
                    Sf += gL.weights[p] * D[q][p] * u[k][p];
                }
                Sf *= a;

                // calculate numerical flux calculation for row
                fLL = numericalFLux[k] * D.LqL[q];         // fluxL(xK) = numericalFlux(k) * lq(xL)
                fLR = numericalFLux[k + 1] * D.LqR[q];     // fluxR(xK) = numericalFlux(k+1) * lq(xR)

                // calculate rhs of linear advection equation (rhsu: cell x row)
                rhsu[k][q] =  ( Sf - (fLR - fLL)) / (J * gL.weights[q]);

            }
        }
    }

    void k(double** val, const storedGLobatto<Np - 1> & gL, const LagrangianPoly<Np>& D)
    {
        NumericalFluxK(val, gL.nodes, D);

        double Sf {};        // (S)^{T} * (au)
        double fLR {};       // f * lq @ xR (within cell k)
        double fLL {};       // f * lq @ xL (within cell k)

        for (unsigned int k = 0; k < K; k++) // iterate over cells
        {
            for (unsigned int q = 0; q < Np; q++) // iterate over lagrangian functions (row of vector)
            {
                Sf = 0; // placeholder value for each row

                // calculate stiffness * f vector (row)
                for (unsigned int p = 0; p < Np; p++) // iterate over nodes
                {
                    Sf += gL.weights[p] * D[q][p] * val[k][p];
                }
                Sf *= a;

                // calculate numerical flux calculation for row
                fLL = numericalFLux[k] * D.LqL[q];         // fluxL(xK) = numericalFlux(k) * lq(xL)
                fLR = numericalFLux[k + 1] * D.LqR[q];     // fluxR(xK) = numericalFlux(k+1) * lq(xR)

                // calculate rhs of linear advection equation (rhsu: cell x row)
                rhsu[k][q] = (Sf - (fLR - fLL)) / (J * gL.weights[q]);
            }
        }
    }

    // getter
    double*& operator[](size_t index) { return rhsu[index]; }
    double*& operator[](size_t index) const { return rhsu[index]; }

//    double copy(size_t i, size_t j)
//    { return rhsu[i][j]; }

};


template <size_t K, size_t Np>
void LSERK4(Solution& u, AdvectionRHS<K, Np>& AdvRHS, const storedGLobatto<Np - 1>& gL, const LagrangianPoly<Np>& LagrangePoly,
         const double& dt)
{
    constexpr static double a[5] = {(double)0,
                                    -(double)567301805773 / 1357537059087,
                                    -(double)2404267990393 / 2016746695238,
                                    -(double)3550918686646 / 2091501179385,
                                    -(double)1275806237668 / 842570457699};

    constexpr static double b[5] = {(double)1432997174477 / 9575080441755,
                                    (double)5161836677717 / 13612068292357,
                                    (double)1720146321549 / 2090206949498,
                                    (double)3134564353537 / 4481467310338,
                                    (double)2277821191437 / 14882151754819};

    constexpr static double c[5] = {(double)0,
                                    (double)1432997174477 / 9575080441755,
                                    (double)2526269341429 / 6820363962896,
                                    (double)2006345519317 / 3224310063776,
                                    (double)2802321613138 / 2924317926251};

    // storage variable
    double pFlat[K * Np] {};
    double* p[K] {};
    double k[K * Np] {};

    // step 0 p(0) = u(n),
    for (size_t i = 0; i < K; i++)
    {
        p[i] = &(pFlat[i * Np]);
        for (size_t j = 0; j < Np; j++)
        {
            p[i][j] = u.copy(i, j);
        }
    }

    // step 1:
    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < Np; j++)
        {
            k[(i * Np) + j] = AdvRHS[i][j] * dt;        // k(1) = dt * advRHS(0)
            p[i][j] += b[0] * k[(i * Np) + j];          // p(1) = p(0) + b(1) * k(1)
        }
    }

    // remaining 4-step iteration
    for (size_t z = 1; z < 5; z++)
    {
        // iterate through cells
        for (size_t i = 0; i < K; i++)
        {
            // iterate through nodes
            for (size_t j = 0; j < Np; j++)
            {
                // k(z)[i][j] = a[z] * k(z-1)[i][j] + dt * uRHS( p[Np * i + j], time + c[z] * dt )

                AdvRHS.k(p, gL, LagrangePoly);  // updating flux calc using intermediate solution p
                k[(i * Np) + j] = (a[z] * k[(i * Np) + j]) + AdvRHS[i][j] * dt; // iterating k term (how to include time?)

                // p(z) = p(z-1) + b[z] * k(z) // iterating p term
                p[i][j] += b[z] * k[(i * Np) + j];
            }
        }
    }

    // update solution (at end of each step)
    for (size_t i = 0; i < K; i++)
    {
        // iterate through nodes
        for (size_t j = 0; j < Np; j++)
        {
            u[i][j] = p[i][j];
        }

        // maintain dirichlet BC (inaccurate if using GL nodes)
        u[0][0] = u.analytical(0, 0);
        u[K - 1][Np - 1] = u.analytical(K - 1, Np - 1);
    }
}


template <size_t K, size_t Np>
void RK4(Solution& u, AdvectionRHS<K, Np>& AdvRHS, const storedGLobatto<Np - 1>& gL, const LagrangianPoly<Np>& LagrangePoly,
            const double& dt)
{
    // storage variables
    double pFlat[K * Np] {};
    double* p[K] {};

    double k1[K * Np] {};
    double k2[K * Np] {};
    double k3[K * Np] {};
    double k4[K * Np] {};

    // step 1:
    for (size_t i = 0; i < K; i++)
    {
        p[i] = &(pFlat[i * Np]);

        for (size_t j = 0; j < Np; j++)
        {
            k1[(i * Np) + j] = AdvRHS[i][j] * dt;        // k(1) = dt * advRHS(u)

            p[i][j] = (u[i][j] + k1[(i * Np) + j] / 2);
        }
    }

    // step 2:
    AdvRHS.k(p, gL, LagrangePoly); // update values for advRHS(u + k1 / 2)
    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < Np; j++)
        {
            k2[(i * Np) + j] = AdvRHS[i][j] * dt;        // k(2) = dt * advRHS(u + k1 / 2)

            p[i][j] = (u[i][j] + k2[(i * Np) + j] / 2);
        }
    }

    // step 3:
    AdvRHS.k(p, gL, LagrangePoly); // update values for advRHS(u + k2 / 2)
    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < Np; j++)
        {
            k3[(i * Np) + j] = AdvRHS[i][j] * dt;        // k(3) = dt * advRHS(u + k2 / 2)

            p[i][j] = (u[i][j] + k3[(i * Np) + j]);
        }
    }

    // step 4:
    AdvRHS.k(p, gL, LagrangePoly); // update values for advRHS(u + k3)
    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < Np; j++)
        {
            k4[(i * Np) + j] = AdvRHS[i][j] * dt;        // k(4) = dt * advRHS(u + k3)
        }
    }

    // update solution (at end of RK4)
    for (size_t i = 0; i < K; i++)
    {
        // iterate through nodes
        for (size_t j = 0; j < Np; j++)
        {
            u[i][j] += ((k1[(i * Np) + j] + k4[(i * Np) + j]) / 6) + ((k2[(i * Np) + j] + k3[(i * Np) + j]) / 3) ;
        }
        // maintain dirichlet BC (to make accurate without solving system must use Gauss-Lobatto Nodes?)
        u[0][0] = u.analytical(0, 0);
        u[K - 1][Np - 1] = u.analytical(K - 1, Np - 1);
    }
}

// export complete graph of solution over domain
template <size_t Np>
void writeL(const size_t& K, const double* gLnodes, const LagrangianPoly<Np>& Lagrange, const Grid& grid, Solution& u)
{
    const int num = 800;  // number of points for plotting per cell (multiply by K for total points)

    // placeholder variables
    double x;
    double r;

    // vertices
    double xL;
    double xR;

    double step; // in this case constant, but leaving general

    // open files to write to
    std::ofstream xFile("Data/xL.txt");
    std::ofstream LFile("Data/L.txt");

    for (size_t i = 0; i < K; i++) // iterate over each cell
    {
        xL = grid.vX[i];
        xR = grid.vX[i + 1];

        x = xL;
        step = (xR - xL) / (num - 1);

        for (size_t j = 0; j < num; j++) // sample at num points
        {
            r = mapAffine(x, xL, xR);
            // removing space at end
            if (j == num - 1 && i == K - 1)
            {
                xFile << x;
                LFile << Lagrange.L(r, gLnodes, u[i]);
            }
            else
            {
                xFile << x << " ";
                LFile << Lagrange.L(r, gLnodes, u[i]) << " ";
            }
            // increment x each time
            x += step;
        }
    }

    LFile << "\n";
    // close files
    xFile.close();
    LFile.close();

}

// naming convention so different iterations aren't overwritten
class Classification
{
private:
    std::string parameters;

public:
    Classification(const unsigned int& N, const double& cfl)
    {
//        parameters = "Cells " + std::to_string(K) + " Order " + std::to_string(N) + " CFL " + std::to_string(cfl);
        parameters = "Order " + std::to_string(N) + " CFL " + std::to_string(cfl);

    }

    [[nodiscard]] std::string read() const
    { return parameters; }

    // writing classification to file
    void write() const
    {
        // writing classification to document
        std::ofstream classifications;
        classifications.open("Data/Iterations.txt", std::ios_base::app);

        if (!classifications.is_open())
        { std::cout << "Error: Could not open file." << std::endl; }
        else { classifications << parameters << " -"; }

        classifications.close();
    }
};


// Solution function | running each iteration -------------------------------------------------------- |
template <size_t OrderN, size_t K>
void solve(const double& cfl, const double& tFinal) // input: order, number of nodes, timestep, domain
{
    // number of nodes in each cell
    const unsigned int Np {OrderN + 1}; // 2 through 11

    // classification (className) of solution
    const Classification className(OrderN, cfl);

//    // extracting Gauss-Legendre nodes and weights
//    const storedGL<OrderN> gL = GaussLegendre<OrderN>();

    // extracting Gauss-Lobatto nodes and weights
    const storedGLobatto<OrderN> gL = GaussLobatto<OrderN>();

    // initializing grid (K x Np)
    const Grid grid(K, Np, gL.nodes);

    const double dx = grid.dx();

    // initializing lagrangian polynomial
    LagrangianPoly<Np> Lagrange(gL.nodes);

    // initializing solution (K x Np)
    Solution u(grid, K, OrderN, Np);

    // initializing advection discretized equation
    AdvectionRHS<K, Np> uRHS(u, grid, gL, Lagrange);

    // initialize timestep using CFL criterion
    const double dt {cfl * dx / a};

    // count number of iterations
    int count {};

//    // write initial condition
//    u.writeSoln(className.read(), count);

//    // plot entire domain
//    writeL(K, gL.nodes, Lagrange, grid, u);

    // iteration
    while (u.Time() < tFinal) {

        u.timestep(dt); // update time
        count++;

        u.analyticalSoln(grid); // update analytical solution

        uRHS.update(u, gL, Lagrange); // update RHS of advection equation
        RK4<K, Np>(u, uRHS, gL, Lagrange, dt); // iterate solution using RK4

        u.residualCalc(); // update residual (should I be using more values on domain instead of just grid pts?)
    }

    // writing solution to file at final time
//    className.write();
////    grid.write(className.read());
////    u.writeSoln(className.read(), count);
//    grid.writeGridSpacing(className.read());
//    u.writeResidual<Np>(className.read(), grid, gL.nodes, Lagrange);

    // plot over entire domain
    writeL(K, gL.nodes, Lagrange, grid, u);
    std::cout << "Final Time: " << u.Time() << std::endl;
}

// main function | defining parameters and running solve functions -----------------------------------|

int main() {

    { // initialize variables:
        const int N{10}; // order of scheme (1 through 10)

        double cfl{0.01}; // Courant number

        const double tFinal {10}; // end time

        // iterate over order:
//        solve<N, 2>(cfl, tFinal);
//        solve<N, 2>(cfl, tFinal);
//        solve<N, 4>(cfl, tFinal);
//        solve<N, 8>(cfl, tFinal);
//        solve<N, 16>(cfl, tFinal);
        solve<N, 1000>(cfl, tFinal);
//        solve<N, 64>(cfl, tFinal);
//        solve<N, 128>(cfl, tFinal);
//        solve<N, 256>(cfl, tFinal);
//        solve<N, 1500>(cfl, tFinal);
    }
        std::cout << "\n" << "Completed" << std::endl;

        return 0;
}
