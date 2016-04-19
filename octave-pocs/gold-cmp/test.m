W = W = full(mmread("1000.mtx"));
D = diag(sum(W));
L = D - W;
tic; [V,l]= eig(L); toc
ac = l(2,2);
fv = V(:,2);
printf("ac = %.16f\n", ac);
fv
res = norm(L*fv - ac*fv, inf) / norm(L, inf);
printf("res = %.16f\n", res);
printf("res.num1 = %.16f, res.num2 = %.16f, res.den=%.16f\n",
       norm(L*fv, inf), norm(ac*fv, inf), norm(L, inf));


