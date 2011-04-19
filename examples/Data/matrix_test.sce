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

fprintfMat("A.txt", A, format="%.18f");
fprintfMat("tA.txt", tA, format="%.18f");
fprintfMat("B.txt", B, format="%.18f");
fprintfMat("tB.txt", tB, format="%.18f");
fprintfMat("C.txt", C, format="%.18f");
fprintfMat("AB.txt", A*B, format="%.18f");
fprintfMat("alpha_AB_beta_C.txt", alpha * A * B + beta * C, format="%.18f");
fprintfMat("alpha_Ax_beta_y.txt", alpha * A * x + beta * y, format="%.18f");

fprintfMat("x.txt", x, format="%.18f");
fprintfMat("y.txt", y, format="%.18f");
fprintfMat("Ax.txt", A * x, format="%.18f");

fprintfMat("cumsum_A_r.txt", cumsum(A, 'r'), format="%.18f");
fprintfMat("cumsum_A_c.txt", cumsum(A, 'c'), format="%.18f");
fprintfMat("cumprod_A_r.txt", cumprod(A, 'r'), format="%.18f");
fprintfMat("cumprod_A_c.txt", cumprod(A, 'c'), format="%.18f");


fprintfMat("sum_A_r.txt", sum(A, 'r')', format="%.18f");
fprintfMat("sum_A_c.txt", sum(A, 'c'), format="%.18f");
fprintfMat("prod_A_r.txt", prod(A, 'r')', format="%.18f");
fprintfMat("prod_A_c.txt", prod(A, 'c'), format="%.18f");

fprintfMat("min_A_r.txt", min(A, 'r')', format="%.18f");
fprintfMat("min_A_c.txt", min(A, 'c'), format="%.18f");
fprintfMat("min_A_star.txt", min(A, '*'), format="%.18f");
fprintfMat("max_A_r.txt", max(A, 'r')', format="%.18f");
fprintfMat("max_A_c.txt", max(A, 'c'), format="%.18f");
fprintfMat("max_A_star.txt", max(A, '*'), format="%.18f");

[ A_sort, ind ] = sort(A, 'c', 'd');
fprintfMat("sort_A_c_d.txt", A_sort, format="%.18f");
fprintfMat("sort_A_c_d_index.txt", ind, format="%.0f");

[ A_sort, ind ] = sort(A, 'c', 'i');
fprintfMat("sort_A_c_i.txt", A_sort, format="%.18f");
fprintfMat("sort_A_c_i_index.txt", ind, format="%.0f");

[ A_sort, ind ] = sort(A, 'r', 'd');
fprintfMat("sort_A_r_d.txt", A_sort, format="%.18f");
fprintfMat("sort_A_r_d_index.txt", ind, format="%.0f");

[ A_sort, ind ] = sort(A, 'r', 'i');
fprintfMat("sort_A_r_i.txt", A_sort, format="%.18f");
fprintfMat("sort_A_r_i_index.txt", ind, format="%.0f");
