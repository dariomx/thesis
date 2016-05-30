using Base.SparseMatrix.CHOLMOD

function calc_fiedler(L)
    n = size(L)[1];
    F = cholfact(L);    
	  ac_upbound = n/(n-1) * minimum(diag(L+1));    
    a = 1;
	  b = (ac_upbound + a)/n;
    v1 = Dense(ones(n) * sqrt(b)); 
    cholmod_updown(1, v1, F);
    mysolve(b) = F \ b;
    lams, vs, _, _, _, _ = eigs(L, nev=1, which=:LM,
                                sigma=0, tol=1e-7,
                                custom_solver=true,
                                mysolve=mysolve);
    ac, fv = lams[1], vs[:,1];
    return ac, fv
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

