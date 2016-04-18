W = W = full(mmread("867e.mtx"));
D = diag(sum(W));
L = D - W;
tic; [V,l]= eig(L); toc
ac = l(2,2);
fv = V(:,2);
printf("ac = %.16f\n", ac);
res = norm(L*fv - ac*fv, inf) / norm(L, inf);
printf("res = %.16f\n", res);



