extern void dsaupd_
(
 int    * ido,
 char   * bmat,
 int    * n,
 char   * which,
 int    * nev,
 double * tol,
 double * resid,
 int    * ncv,
 double * V,
 int    * ldv,
 int    * iparam,
 int    * ipntr,
 double * workd,
 double * workl,
 int    * lworkl,
 int    * info
);

extern void dseupd_
(
 bool   * rvec,
 char   * howmny,
 bool   * select,
 double * d,
 double * Z,
 int    * ldz,
 double * sigma,
 char   * bmat, 
 int    * n,
 char   * which,
 int    * nev,
 double * tol,
 double * resid,
 int    * ncv,
 double * V,
 int    * ldv,
 int    * iparam,
 int    * ipntr,
 double * workd,
 double * workl,
 int    * lworkl,
 int    * info 
);
