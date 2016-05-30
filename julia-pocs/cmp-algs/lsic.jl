using Base.SparseMatrix.CHOLMOD

function calc_fiedler(L)
    a = 0
	  La = L + a*I;
    F = cholfact(La);
    mysolve(b) = F \ b;
    lams, vs, _, _, _, _ = eigs(L, nev=2, which=:LM,
                                sigma=0, tol=1e-7,
                                custom_solver=true,
                                mysolve=mysolve);
    ac, fv = lams[2], vs[:,2];
    return ac + a, fv
end    

function test_fiedler(fn)
    L = CHOLMOD.Sparse(fn);
    @time ac, fv = calc_fiedler(L);
    mynorm(x) = norm(x, Inf);
    res = mynorm(L*fv - ac*fv);
    @printf("ac=%.14f, res=%.3E\n", ac, res);
end

# main
fn = ARGS[1];
test_fiedler(fn);
test_fiedler(fn);
test_fiedler(fn);
