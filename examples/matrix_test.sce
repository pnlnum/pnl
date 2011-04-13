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
