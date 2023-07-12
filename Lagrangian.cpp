// Created by Aidan Christoffersson on 2023-06-13.

#include <iostream>
#include <cmath>
#include <fstream>

// define Lagrangian interpolant polynomial
template <size_t Np>
class LagrangianPoly
{
private:
    // constant denominator of L
    double denom[Np] {};

    // storing derivative of L at all nodes
    const std::unique_ptr<double[]> DFlattened;
    const std::unique_ptr<double*[]> D;

public:
    // boundary values for basis functions | as nodes const location -> constant functions, no?
    double LqL[Np] {};
    double LqR[Np] {};

private:
    double dL(const double* gLnodes, size_t& q, size_t& p) // computes the derivative of L[q] at x[p]
    {
        double dL {0};

        for (unsigned int l = 0; l < Np; l++) // iterate over sum (result of chain rule)
            { if (l != q) // q is the index of the Lagrange poly
                { double k = pow((gLnodes[q] - gLnodes[l]), -1);
                    for (unsigned int i = 0; i < Np; i++)
                    { if ((i != q) && (i != l)) // avoid division by zero
                        {
                            k *= (gLnodes[p] - gLnodes[i]) / (gLnodes[q] - gLnodes[i]);
                        }
                    }
                    // add value to sum
                    dL += k;
                }
            }
        return dL;
    }

    void LqUpdate(const double* gLnodes)
    {
        // boundary of cell domain
        double rR {1};
        double rL {-1} ;

        for (size_t q = 0; q < Np; q++) // iterate through each Lagrangian function (coln)
        {
            LqL[q] = Lq(rL, gLnodes, q);
            LqR[q] = Lq(rR, gLnodes, q);
        }
    };

public:

    double Lq(double r, const double* gLnodes, size_t& q) const // returns value of Lagrangian polynomial q at point r
    {
        double Lq {1};
        for (size_t p = 0; p < Np; p++) // iterating over each r node to calculate each L[q]
        {
            if (p != q)
            {
                Lq *= (r - gLnodes[p]); // product of all numerator terms
            }
        }
        return (Lq / denom[q]);
    }

    //constructor for Lagrangian interpolating polynomial
    explicit LagrangianPoly(const double* gLnodes)
    : DFlattened(std::make_unique<double[]>(Np * Np)), D(std::make_unique<double*[]>(Np))
    {

        // initializing denominator
        for (unsigned int q = 0; q < Np; q++) // array is the set of Li from 0 -> Np - 1
        {
            denom[q] = 1; // holding variable for products

            for (unsigned int p = 0; p < Np; p++) // compute each Li
            {
                if (p != q)
                { denom[q] *= (gLnodes[q] - gLnodes[p]); }
            }
        }

        // calculating derivative matrix
        for (size_t q = 0; q < Np; q++) // q is the Lagrangian index
        {
            D[q] = &(DFlattened[q * Np]);

            for (size_t p = 0; p < Np; p++) // p is the node index
            {
                D[q][p] = dL(gLnodes, q, p);
            }
        }

        // calculating end nodes for lagrangian basis functions
        LqUpdate(gLnodes);

//        // printing D matrix
//        for (unsigned int q = 0; q < Np; q++) // q is the Lagrangian index
//        {
//            for (unsigned int p = 0; p < Np; p++) // p is the node index
//            {
//                std::cout << D[q][p] << "  ";
//            }
//            std::cout << std::endl;
//        }
    }

    // returns value of Lagrangian polynomial at r | used for extrapolation to find value at endpoints
    double L(double r, const double* gLnodes, const double* uCell) const
    { // use affine mapping in function implementation to project from x to r, uCell is ptr to array of cell nodes
        double soln {};
        for (size_t q = 0; q < Np; q++) // summing up each Lagrangian Polynomial multiplied by weight. (soln node)
        {
            soln += Lq(r, gLnodes, q) * uCell[q];
        }
        return soln;
    }

    // getter
    double*& operator[](size_t index) { return D[index]; }
    double*& operator[](size_t index) const { return D[index]; }

    // write Lagrangian Polynomials to file for testing
    void testLq(const double* gLnodes)
    {
        const size_t num = 500;
        double r;

        const double step {(double)2 / ((double)num - 1)};

        // open files to write to
        std::ofstream rL("Data/rL.txt");
        std::ofstream Lfile("Data/L.txt");

        for (size_t q = 0; q < Np; q++) // iterate through each Lagrangian poly
        {
            r = -1;

            for (int i = 0; i < num; i++) // iterate through domain for plotting
            {
                // removing space at end
                if (i == num - 1)
                {
                    rL << r;
                    Lfile << Lq(r, gLnodes, q);
                }
                else
                {
                    rL << r << " ";
                    Lfile << Lq(r, gLnodes, q) << " ";
                }

               // update r
               r += step;
            }

            rL << std::endl;
            Lfile << std::endl;

        }

        // close files
        rL.close();
        Lfile.close();

    }

};