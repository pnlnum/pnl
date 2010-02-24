#include "blaswrap.h"
#include "f2c.h"

/* Subroutine */ int cspsv_(char *uplo, integer *n, integer *nrhs, complex *
	ap, integer *ipiv, complex *b, integer *ldb, integer *info)
{
/*  -- LAPACK driver routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    CSPSV computes the solution to a complex system of linear equations   
       A * X = B,   
    where A is an N-by-N symmetric matrix stored in packed format and X   
    and B are N-by-NRHS matrices.   

    The diagonal pivoting method is used to factor A as   
       A = U * D * U**T,  if UPLO = 'U', or   
       A = L * D * L**T,  if UPLO = 'L',   
    where U (or L) is a product of permutation and unit upper (lower)   
    triangular matrices, D is symmetric and block diagonal with 1-by-1   
    and 2-by-2 diagonal blocks.  The factored form of A is then used to   
    solve the system of equations A * X = B.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    AP      (input/output) COMPLEX array, dimension (N*(N+1)/2)   
            On entry, the upper or lower triangle of the symmetric matrix   
            A, packed columnwise in a linear array.  The j-th column of A   
            is stored in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   
            See below for further details.   

            On exit, the block diagonal matrix D and the multipliers used   
            to obtain the factor U or L from the factorization   
            A = U*D*U**T or A = L*D*L**T as computed by CSPTRF, stored as   
            a packed triangular matrix in the same storage format as A.   

    IPIV    (output) INTEGER array, dimension (N)   
            Details of the interchanges and the block structure of D, as   
            determined by CSPTRF.  If IPIV(k) > 0, then rows and columns   
            k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1   
            diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,   
            then rows and columns k-1 and -IPIV(k) were interchanged and   
            D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and   
            IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and   
            -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2   
            diagonal block.   

    B       (input/output) COMPLEX array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization   
                  has been completed, but the block diagonal matrix D is   
                  exactly singular, so the solution could not be   
                  computed.   

    Further Details   
    ===============   

    The packed storage scheme is illustrated by the following example   
    when N = 4, UPLO = 'U':   

    Two-dimensional storage of the symmetric matrix A:   

       a11 a12 a13 a14   
           a22 a23 a24   
               a33 a34     (aij = aji)   
                   a44   

    Packed storage of the upper triangle of A:   

    AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* System generated locals */
    integer b_dim1, b_offset, i__1;
    /* Local variables */
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int xerbla_(char *, integer *), csptrf_(
	    char *, integer *, complex *, integer *, integer *), 
	    csptrs_(char *, integer *, integer *, complex *, integer *, 
	    complex *, integer *, integer *);

    --ap;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CSPSV ", &i__1);
	return 0;
    }

/*     Compute the factorization A = U*D*U' or A = L*D*L'. */

    csptrf_(uplo, n, &ap[1], &ipiv[1], info);
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

	csptrs_(uplo, n, nrhs, &ap[1], &ipiv[1], &b[b_offset], ldb, info);

    }
    return 0;

/*     End of CSPSV */

} /* cspsv_ */

