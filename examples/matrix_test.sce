m = 7;
n= 5;

A = rand (m, n);
tA = A';
B = rand (n, m);
tB = B'
C = rand (m,  m);
x = rand (n,1);
y = rand (m,1);
alpha = 1.3;
beta = 0.5;

fprintfMat("Data/A.txt", A, format="%.18f");
fprintfMat("Data/tA.txt", tA, format="%.18f");
fprintfMat("Data/B.txt", B, format="%.18f");
fprintfMat("Data/tB.txt", tB, format="%.18f");
fprintfMat("Data/C.txt", C, format="%.18f");
fprintfMat("Data/AB.txt", A*B, format="%.18f");
fprintfMat("Data/alpha_AB_beta_C.txt", alpha * A * B + beta * C, format="%.18f");
fprintfMat("Data/alpha_Ax_beta_y.txt", alpha * A * x + beta * y, format="%.18f");

fprintfMat("Data/x.txt", x, format="%.18f");
fprintfMat("Data/y.txt", y, format="%.18f");
fprintfMat("Data/Ax.txt", A * x, format="%.18f");

fprintfMat("Data/cumsum_A_r.txt", cumsum(A, 'r'), format="%.18f");
fprintfMat("Data/cumsum_A_c.txt", cumsum(A, 'c'), format="%.18f");
fprintfMat("Data/cumprod_A_r.txt", cumprod(A, 'r'), format="%.18f");
fprintfMat("Data/cumprod_A_c.txt", cumprod(A, 'c'), format="%.18f");


fprintfMat("Data/sum_A_r.txt", sum(A, 'r')', format="%.18f");
fprintfMat("Data/sum_A_c.txt", sum(A, 'c'), format="%.18f");
fprintfMat("Data/prod_A_r.txt", prod(A, 'r')', format="%.18f");
fprintfMat("Data/prod_A_c.txt", prod(A, 'c'), format="%.18f");

fprintfMat("Data/min_A_r.txt", min(A, 'r')', format="%.18f");
fprintfMat("Data/min_A_c.txt", min(A, 'c'), format="%.18f");
fprintfMat("Data/max_A_r.txt", max(A, 'r')', format="%.18f");
fprintfMat("Data/max_A_c.txt", max(A, 'c'), format="%.18f");

