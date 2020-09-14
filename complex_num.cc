//Author: Deng Pan, Scarlett Gong

#include <stdio.h>
#include <complex>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>

using namespace std; 
using Complex = complex<double>;

Complex pz(Complex guess){
    //define the polynomial p(z) = z^3-1
    Complex z = pow(guess,3) - 1.;
    return z;
}

Complex pzDiff(Complex guess){
    //define the derivative of p(z) = z^3-1 
    Complex z = 2. *pow(guess,2) ;
    return z;
}
Complex iterNewton(Complex guess, Complex (*poly)(Complex), Complex (*polyDiff)(Complex)){
    //one iteration Newton's method
    Complex z = guess - poly(guess)/polyDiff(guess);
    return z;
}
int newtonFractal(Complex initial_guess, Complex (*poly)(Complex), Complex (*polyDiff)(Complex), Complex *roots , double tol, int max_iterations){
    //Given a initial_guess, a polynomial function and its derivative, and the roots we are expecting it to converge, return the Int that represents the converged roots
    //return 0: does not converge, return other positive integer i : the ith roots in roots list 
    int n_roots = sizeof(roots);
    Complex z = iterNewton(initial_guess, poly, polyDiff);
    //for each iteration if some root is close to result (less than tolerance distance) return that root index otherwise return 0 after maximum iterations.
    for (int i = 0; i<max_iterations;i++){
        for (int j = 0; j<n_roots; j++){
            if (abs(roots[j]- z ) < tol) {
                return j+1;
            }
        }
        z = iterNewton(z, poly, polyDiff);
    }
    return 0;
}
int main(){
    clock_t t = clock();
    using namespace std;

    // Half length of the subject of complex plane we will
    // represent
    double L = 4;
    // The number of grid points along one dimension
    int N = 1000;
    // The grid spacing
    double dh = 2*L/( (double) (N-1));

    // declare and initialize the complex number at the lower left corner
    complex<double> z_corner (-L , -L);
    // declare an 1D array, which we will use index manipulation to
    // represent a 2D array
    complex<double>  * initial_guess = new complex<double> [N*N];

    Complex roots[3] = {Complex(1,0), Complex(-0.5, sqrt(3)/2), Complex(-0.5 , - sqrt(3)/2)};
    int max_iters = 1000;
    double tol = 1e-6;

    // Imagine a 2D NxN array, we take the second row and put it behind
    // the first row, and take the third and put it behind the second row,
    // so on and so forth
    // In the end, we have 1D array of length NxN
    // To get the (j,i) element in the original 2D array,
    // we will need to multiply the row index by N, and add the column index,
    // i.e. index_{1D} = j*N + i.
    ofstream data_file("data_file.txt");
    for (int j=0;j<N;j++) {
        for (int i=0;i<N;i++) {
            initial_guess [i+j*N] = z_corner + complex<double>(dh*i, dh*j);
            data_file << real(initial_guess[i+j*N]) << " " <<imag(initial_guess[i+j*N]) << " "<< newtonFractal(initial_guess[i+j*N],pz,pzDiff, roots, tol, max_iters) << endl;
        }
    }

    delete [] initial_guess;
    t = clock() - t;
    printf("We use %f seconds for calculation of size N = %d  \n", (float)t/CLOCKS_PER_SEC, N);
}
