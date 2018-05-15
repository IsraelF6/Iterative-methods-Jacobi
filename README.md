# Iterative-methods-Jacobi


A subroutine that solves a triadiagonal matrix system, Ax = b, using the Jacobi method. 
Your subroutine should take as input the matrix coefficients along the 3 bands: 
a[i,i−1] where i=2, ..., n; a[i,i] where i=1, ..., n, and a[i,i+1] where i=1, ..., n−1. 
Also b, x(0), and the tolerance eps, are inputs.
Your subroutine should only store the nonzero entries (the 3 bands of A)
of A. Compare your new Jacobi method subroutine (where only the 3
bands of A are stored) to the Gaussian Elimination approach (where the
whole matrix A, including zero values, is stored) in regards to efficiency
and accuracy when the n × n matrix A and the n × 1 vector b are defined
as,

a[i,j] =			{ 4 if i=j,  −1 if j=i+1,  −1 if j=i−1,  0 otherwise

b[i] =			{ 1 if i<n/2,  0 otherwise }

eps = 10^-10 and x(0)=0.
