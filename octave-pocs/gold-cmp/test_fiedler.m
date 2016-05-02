#!/usr/bin/octave -qf

warning('off',
        'Octave:possible-matlab-short-circuit-operator');
args = argv();
fn = args{1};
printf("reading weight matrix from %s\n", fn)
W = full(mmread(fn));
D = diag(sum(W));
L = D - W;

printf("calculating alg-conn and fiedler vector ...\n");
tic; [V,l]= eig(L); toc
ac = l(2,2);
fv = V(:,2);

printf("ac = %.16f\n", ac);
#fv
res = norm(L*fv - ac*fv, inf) / norm(L, inf);
printf("res = %.16f\n", res);
printf("res.num1 = %.16f\nres.num2 = %.16f\nres.den = %.16f\n",
       norm(L*fv, inf), norm(ac*fv, inf), norm(L, inf));
