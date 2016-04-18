W = W = full(mmread("5000.mtx"));
D = diag(sum(W));
L = D - W;
numEvals = 2;
target = 'SA';
primme_start = tic;
[primme_V,primme_D,norms,primmeout] = primme_eigs(L, numEvals, target, struct(), eigsMethod=0);
primme_time_elapsed = toc(primme_start);
ac = diag(primme_D)(2:2);
fv = primme_V(:,2);
relres = norm(L*fv - ac*fv, inf) / norm(L,inf);
printf ("calc took: %.16f\n", primme_time_elapsed);
printf ("alg conn: %.16f\n", ac);
printf("relres c: %.16f\n", relres);




