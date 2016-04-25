printf ("reading matrix from file ..."); fflush(stdout);
#L = full(mmread("domain/867.mtx"));
W = full(mmread("mine/6000.mtx"));
printf ("calculating laplacian ..."); fflush(stdout);
D = diag(sum(W));
L = D - W;
numEvals = 2;
target = 'SA';
method = 0;
printf ("calling eigen solver ..."); fflush(stdout);
primme_start = tic;
[primme_V,primme_D,norms,primmeout] = primme_eigs(L, numEvals, target, struct("precontition",1), eigsMethod=method);
primme_time_elapsed = toc(primme_start);
printf ("calculating residuals ..."); fflush(stdout);
ac = diag(primme_D)(2:2);
fv = primme_V(:,2);
relres = norm(L*fv - ac*fv, 2) / norm(L, 2);
printf ("calc took: %.16f\n", primme_time_elapsed);
printf ("alg conn: %.16f\n", ac);
printf("relres c: %.16f\n", relres);




