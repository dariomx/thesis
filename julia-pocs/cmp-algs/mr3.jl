using Base.SparseMatrix.CHOLMOD

function calc_fiedler(L)
    lams, vs = eig(L)
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

