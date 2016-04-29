!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Computes the Fiedler vector using the TraceMin algorithm                !
!                                                                            !
!    Copyright (C) 2010 Murat Manguoglu, Faisal Saied, Eric Cox, Ahmed Sameh,!
!    Purdue University                                                       !
!    Contact: murat.manguoglu@gmail.com                                      !
!                                                                            !
!    This file is a part of TraceMin-Fiedler                                 !
!    TraceMin-Fiedler is free software: you can redistribute it and/or mofy!
!    it under the terms of the GNU General Public License as published by    !
!    the Free Software Foundation, either version 3 of the License, or       !
!    (at your option) any later version.                                     !
!                                                                            !
!    This program is distributed in the hope that it will be useful,         !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of          !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           !
!    GNU General Public License for more details.                            !
!                                                                            !
!    You should have received a copy of the GNU General Public License       !
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.   !
!								             !
!    By downloading the software TraceMin-Fiedler I agree to give reference  !
!    to the authors of the software in any publications resulting from its   !
!    use.                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine FIEDLER(neqns, nzmax, ia, ja, a,fiedler_vec,max_trmin_iters,max_pspike_iters,tracemin_tol,pspike_tol, abs_res_tol, lam_tol, num_procs, rank )
use mpi 
use omp_lib


implicit none 
!include 'omp_lib.h'


double precision  :: a(nzmax) 
integer :: ia(neqns+1) , ja(nzmax)
double precision  :: fiedler_vec(neqns) 
double precision,   dimension(:), allocatable ::Av, x, bs 
double precision ::  aa 
integer,           dimension(:), allocatable ::   Ai, Aj, rcounts, displs  

double precision   ::  c1, c2
character(len=100) :: name,name_rhs

INTEGER*8 pt(64)
INTEGER   mtype
INTEGER   maxfct, mnum, phase, msglvl, nzmax
INTEGER   iparm(64),  spikeparm(100) 
double precision ddum(23), t0 , t1 , bnorm, tol , norma, sqrt  
integer  neqns, i, j, k, error, job, bandwidth ,  n2 , m , n , nnz


INTEGER status2(MPI_STATUS_SIZE, 4)
integer, dimension(MPI_STATUS_SIZE) :: status
integer :: rank, num_procs,code, req(4), ierr

      character rep*10
      character field*7
      character symm*19
      character tmp1*1024
      character tmp2*2
      integer,dimension(:),allocatable   :: IDummy
      double complex, dimension(:),allocatable   :: CDummy

logical decay
logical prec 
logical knownvector 
!===========================================================
! Tracemin variables
!
integer max_trmin_iters, trmin_iters, trmin_iters2, trmin_flag
integer max_pspike_iters
double precision,  dimension(:,:), allocatable :: W, Y, AY, RITZ, WORK, Yconv, b 
double precision,  dimension(:), allocatable :: binout, rjj

integer  ndesired
double precision myrand, seed
integer iseed
integer nrhs, nloc, kk  

integer qr_Lwork, qr_info, lu_LWORK, lu_info, LDschur
integer, dimension(:), allocatable :: qr_jpvt, lu_IPIV
double precision, dimension(:), allocatable :: qr_work, qr_tau
double precision ddot, dnrm2
integer ldz2, info2, index
character*1  jobz, uplo
double precision, allocatable, dimension(:) :: W2, work2, lu_WORK, H_loc, H_glob
double precision, allocatable, dimension(:,:) :: Z2, Schur
double precision  Tstart, T_preprocess1, T_preprocess2, T_solve, T_other
double precision  T_ritz, T_conv, T_orthog, T_defl, Tstart2, T_trmin, Tstart3
double precision  T_trmin_1, T_trmin_2, T_trmin_3, T_trmin_4
double precision  T_ritz_1, T_ritz_2, T_ritz_3, T_ritz_4, T_conv_1, T_conv_2
double precision  T_cgmatvec , T_cg_allred 
integer nconv, ldBmat
double precision  tracemin_tol, pspike_tol, evalue, xjjnorm, abs_res_tol, lam_tol
integer iprob, nx, part1, part1_LAST, jj, j2, Jindex 
double precision,  dimension(:), allocatable :: xtmp, xtmp_loc, lam, resid_norm
double precision,  dimension(:), allocatable :: Schur_loc_buf, Schur_glob_buf, B_loc, B_glob
double precision,  dimension(:,:), allocatable :: Bmat
integer, dimension(:), allocatable :: Jvector
character*1  trans1, trans2
double precision alpha, beta, min_res_norm
integer ncalibrate
double precision T_calibrate_1, T_calibrate_2
integer cgmaxit, cgiter, inner_solver_flag
double precision cgtol, shift
logical debug 
double precision  min_lam_curr, min_lam_prev, min_lam_rel_change
!
!
!===========================================================


!tracemin_tol = 1e-5
!pspike_tol = 1e-4

debug = .true. 

if (rank .eq. 0) then

   ndesired = 2
   inner_solver_flag =2  
   prec = .true. 
   knownvector = .true. 
!   knownvector = .false. 

endif

call MPI_BCAST(ndesired,           1, MPI_INTEGER,0,MPI_COMM_WORLD,error)
call MPI_BCAST(max_trmin_iters,    1, MPI_INTEGER,0,MPI_COMM_WORLD,error)
call MPI_BCAST(max_pspike_iters,   1, MPI_INTEGER,0,MPI_COMM_WORLD,error)
call MPI_BCAST(inner_solver_flag,  1, MPI_INTEGER,0,MPI_COMM_WORLD,error)
call MPI_BCAST(tracemin_tol,       1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,error)
call MPI_BCAST(pspike_tol,         1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,error)
call MPI_BCAST(prec,               1, MPI_LOGICAL,0,MPI_COMM_WORLD,error)
call MPI_BCAST(knownvector,               1, MPI_LOGICAL,0,MPI_COMM_WORLD,error)



nrhs = 3*ndesired 
nconv = 0

min_lam_rel_change = 1000
min_lam_prev = 1000000
!
!===========================================================
   part1 = neqns / num_procs
   part1_LAST = neqns - (num_procs -1 )*part1

   nloc = part1
   if (rank .eq. num_procs-1) nloc = part1_LAST

   allocate(b(nloc,nrhs))
   allocate(Y(nloc,nrhs))
   allocate(W(nloc,nrhs))
   allocate(AY(nloc,nrhs))
   allocate(RITZ(nloc,nrhs));      allocate(WORK(nloc,nrhs))
   allocate(Yconv(nloc,nrhs))


   allocate(xtmp(nrhs), xtmp_loc(nrhs), lam(nrhs), resid_norm(nrhs), Jvector(nrhs))

   allocate(Schur(nrhs,nrhs), Schur_loc_buf(nrhs*nrhs), Schur_glob_buf(nrhs*nrhs))

   allocate (lu_WORK(nrhs))
   LDschur = nrhs
   lu_LWORK = nrhs
   allocate(lu_IPIV(nrhs))



   ldz2=nrhs
   allocate (W2(nrhs), Z2(nrhs,nrhs), work2(3*nrhs))
   allocate(H_glob(nrhs*nrhs), H_loc(nrhs*nrhs))
   ldBmat = nrhs
   allocate(B_loc(nrhs*nrhs), B_glob(nrhs*nrhs), Bmat(nrhs,nrhs))


!
!==================================================================================
!
!  Initialize Y with random values
!

   if (debug .and. (rank .eq. 0)) print *, ' calling random number generator '

   iseed = 107
   seed = dfloat(iseed)/1024.0d0

   do j = 1, rank*part1*nrhs
      seed = myrand(seed)      
   enddo
      
   do i = 1, nloc
      do j = 1, nrhs
         seed = myrand(seed)
         Y(i,j) = seed
      enddo
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!Set the known vector of 1s !
!!!!!Does it really work ??     ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (knownvector) then 
   nconv = 1 
   if (rank .eq. 0) then
      print *, ' knownvector = true '
      print *, ' Starting witn nconv = ', nconv
      print *, ' nrhs = ', nrhs
   endif
   do i = 1,nloc 
      Yconv(i,1) = 1.0d0 / sqrt(float(neqns)) 
!!!???!!!      Y(i,1)  = sqrt(float(neqns)) 
   end do
end if 

call MPI_BCAST(nconv , 1, MPI_INTEGER,0,MPI_COMM_WORLD,error)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------
! deflate initial guess before orthonalizing
if (nconv .ge. 1) then
   Tstart2 = MPI_Wtime()
   call apply_deflation(nloc,nrhs,nconv,Y,Yconv,B_loc,B_glob, Bmat, ldBmat)
   T_defl = T_defl + MPI_Wtime() - Tstart2
endif  !  if (nconv .ge. 1) then (deflate)

T_other = T_other + MPI_Wtime() - Tstart
!-------------------


   if (debug .and. (rank .eq. 0)) print *, ' calling par_orthogonal_basis '

   call par_orthogonal_basis(nloc, nrhs, Y, xtmp, xtmp_loc)

   if (debug .and. (rank .eq. 0)) print *, ' returned from par_orthogonal_basis '
!
!==================================================================================
!==================================================================================
!==================================================================================
!   PSPIKE preprocessing calls (job =0, job = 1)
!
tol = pspike_tol

!
!   The following is for the symmetric case only
!
   Tstart = MPI_Wtime()
      job = 1;
      if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
      if (debug .and. (rank .eq. 0)) print *, ' calling MATVEC SETUP: job = 1 '
      call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,error)
      call MATVEC(job, neqns, nzmax, ia, ja, a, b,nloc, nrhs, num_procs, rank )

      if (rank.eq.0) print *, ' Completed MATVEC SETUP (job = 1) :', rank   

      norma = b(1,1)  
 
      call MPI_BCAST(norma , 1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   T_preprocess2 = MPI_Wtime() - Tstart


T_solve = 0.0d0
T_other = 0.0d0
T_trmin = 0.0d0; T_ritz = 0.0d0; T_conv = 0.0d0; T_orthog = 0.0d0; T_defl = 0.0d0
T_trmin_1 = 0.0d0; T_trmin_2 = 0.0d0; T_trmin_3 = 0.0d0; T_trmin_4 = 0.0d0
T_ritz_1 = 0.0d0; T_ritz_2 = 0.0d0; T_ritz_3 = 0.0d0; T_ritz_4 = 0.0d0
T_conv_1 = 0.0d0; T_conv_2 = 0.0d0
T_cgmatvec = 0.0d0; T_cg_allred = 0.0d0 
if (debug) then 
if (rank .eq. 0)  then
   print *, ' performing ', max_trmin_iters, &
         ' iterations of Trace mininmization'
   print *, ' tracemin_tol = ', tracemin_tol
   print *, ' pcg_tol = ', pspike_tol
   print *, ' Number of desired eigenvalues = ', ndesired
   print *, ' N = ', neqns, ',   bandwidth = ', bandwidth
   print *, ' nloc = ', nloc, ',  nrhs = ', nrhs
   print *, ' norma =   ', norma
endif
end if 
!
!==================================================================================
!

call MPI_BARRIER(MPI_COMM_WORLD,error)

do trmin_iters = 1, max_trmin_iters
   !-------------------
   !
   Tstart = MPI_Wtime()
   if (inner_solver_flag .eq. 1) then ! call spikepardiso

   elseif (inner_solver_flag .eq. 2) then ! call conjugate gradient (modified)
      !===
      cgmaxit = max_pspike_iters
      cgtol = pspike_tol
      !prec = .true.  
      !!!!if (trmin_iters .eq. 1) then 
      if (trmin_iters .gt. 0) then !! this choice says always solve shifted system
         call pcgsingular(prec, neqns, nzmax, ia, ja, a, Y, W, nloc, nrhs, bandwidth, cgtol, &
                 spikeparm, num_procs, rank, xtmp, xtmp_loc, &
                 b, AY, RITZ, WORK, cgiter, cgmaxit, T_cgmatvec, T_cg_allred)
         if (debug .and.(rank .eq. 0 )) &
         write(*,612)  trmin_iters, cgiter
      else 
         call pcg(prec, neqns, nzmax, ia, ja, a, Y, W, nloc, nrhs, bandwidth, cgtol, &
                 spikeparm, num_procs, rank, xtmp, xtmp_loc, &
                 b, AY, RITZ, WORK, cgiter, cgmaxit, T_cgmatvec, T_cg_allred)   
         if (debug .and. (rank .eq. 0 )) &
           write(*,612)  trmin_iters, cgiter
      end if 
   elseif (inner_solver_flag .eq. 3) then ! call conjugate gradient (modified)
      !===
      cgmaxit = max_pspike_iters
      cgtol = pspike_tol
      call pcr(prec, neqns, nzmax, ia, ja, a, Y, W, nloc, nrhs, bandwidth, cgtol, &
                 spikeparm, num_procs, rank, xtmp, xtmp_loc, &
                 b, AY, RITZ, WORK, cgiter, cgmaxit, T_cgmatvec, T_cg_allred)
      if (debug .and. (rank .eq. 0 )) &
           write(*,612)  trmin_iters, cgiter
   else

      if (debug) print *, ' inner_solver_flag must be 1 or 2: ', inner_solver_flag
   endif

611 format(' TraceMIN iter = ', i6, ',   PSPIKE iters performed = ', i5)
612 format(' TraceMIN iter = ', i6, ',   CG iters performed = ', i5)


!do j=1,neqns 
!   write(14,*) Y(j,2) 
!   write(15,*) W(j,2) 
!end do 
!call MPI_BARRIER(MPI_COMM_WORLD,error) 
!stop 

   T_solve = T_solve + MPI_Wtime() - Tstart
   Tstart = MPI_Wtime()
   !-------------------
!!!   Y = W
   !-------------------
   ! Schur complement for Tracemin
     Tstart2 = MPI_Wtime()
!      trmin_flag=0
      trmin_flag=1

      if (trmin_flag .eq. 1) then
!        Schur = Y'*W;         
         Tstart3 = MPI_Wtime()
         k = 1
         do i = 1, nrhs
            do j = 1, i
               ! changed n to nloc
               Schur(i,j)=ddot(nloc,Y(1,i),1,W(1,j),1)
               Schur_loc_buf(k) = Schur(i,j)
               k = k+1
            enddo
         enddo
         T_trmin_1 = T_trmin_1 + MPI_Wtime() - Tstart3

         Tstart3 = MPI_Wtime()
         call MPI_Allreduce(Schur_loc_buf, Schur_glob_buf, k-1, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)
         T_trmin_2 = T_trmin_2 + MPI_Wtime() - Tstart3

         k = 1
         do i = 1, nrhs
            do j = 1, i
               Schur(i,j) = Schur_glob_buf(k)
               k = k+1
               if (i .ne. j) Schur(j,i) = Schur(i,j)
            enddo
         enddo

!        Schur = inv(Schur);
         Tstart3 = MPI_Wtime()
         lu_info=0
         call dgetrf(nrhs, nrhs, Schur, LDschur, lu_IPIV, lu_INFO)
         if (debug .and. (lu_info .ne.0)) print *, 'dgetrf: info = ',lu_info
         lu_info=0
         call dgetri(nrhs, Schur, LDschur, lu_IPIV, lu_WORK, lu_LWORK, lu_INFO)
         if (debug .and. (lu_info .ne.0)) print *, 'dgetri: info = ',lu_info         
         T_trmin_3 = T_trmin_3 + MPI_Wtime() - Tstart3
            
!        Y = W*Schur;         
         Tstart3 = MPI_Wtime()

!!!!    C := alpha*op( A )*op( B ) + beta*C,
!!!!    Y = 1 * W * Schur + 0 * Y
!!!!    DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

    trans1='N';  trans2='N'
    alpha = 1.0d0;   beta = 0.0d0
    call dgemm(trans1,trans2,nloc,nrhs,nrhs,alpha,W,nloc,Schur,LDschur,beta,Y,nloc)


!         do i = 1, nloc 
!            do j = 1, nrhs
!               Y(i,j)=0.0d0
!               do k = 1, nrhs
!                  Y(i,j)=Y(i,j)+W(i,k)*Schur(k,j)
!               enddo
!            enddo
!         enddo
         T_trmin_4 = T_trmin_4 + MPI_Wtime() - Tstart3
      else if (trmin_flag .eq. 0) then
         Y = W
      else
         if (debug .and. (rank .eq. 0)) print *, ' wrong trmin_flag: ', trmin_flag
      endif

   T_trmin = T_trmin + MPI_Wtime() - Tstart2

   !-------------------
   ! deflate
   ! (moved here for test)
   if (nconv .ge. 1) then
      Tstart2 = MPI_Wtime()
      call apply_deflation(nloc,nrhs,nconv,Y,Yconv,B_loc,B_glob, Bmat, ldBmat)
      T_defl = T_defl + MPI_Wtime() - Tstart2
   endif  !  if (nconv .ge. 1) then (deflate)

   T_other = T_other + MPI_Wtime() - Tstart
   !-------------------



   !-------------------
   Tstart2 = MPI_Wtime()
   call par_orthogonal_basis(nloc, nrhs, Y, xtmp, xtmp_loc)
   T_orthog = T_orthog + MPI_Wtime() - Tstart2
   !-------------------
! --- Rayleigh-Ritz acceleration -------------
!

!   if (trmin_iters .le. 4) goto 9090
!   goto 9090

   Tstart2 = MPI_Wtime()
!....... AY = A*Y;
   Tstart3 = MPI_Wtime()
   do j = 1, nrhs
      do i = 1, nloc; AY(i,j) = Y(i,j); enddo
   enddo

   job = 55; !! JOB = 2 for solve 
   if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
   call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MATVEC(job, neqns, nzmax, ia, ja, a, AY,nloc, nrhs, num_procs, rank )
   T_ritz_1 = T_ritz_1 + MPI_Wtime() - Tstart3
!....... H = Y'*A*Y = Y'*AY;
   Tstart3 = MPI_Wtime()
   kk = 1
   do i = 1, nrhs
      do j = i, nrhs
         index = i + (j-1)*j/2
         H_loc(index) = ddot(nloc,Y(1,i),1,AY(1,j),1)
         kk = kk + 1
      enddo
   enddo

   call MPI_Allreduce(H_loc, H_glob, kk-1, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)
   T_ritz_2 = T_ritz_2 + MPI_Wtime() - Tstart3
!....... Solve eigenvalue problem         
   Tstart3 = MPI_Wtime()
   jobz = 'V';      uplo = 'U';      info2 = 0
   call dspev(jobz, uplo, nrhs, H_glob, W2, Z2, ldz2, work2, info2)
   T_ritz_3 = T_ritz_3 + MPI_Wtime() - Tstart3
!....... Form Ritz vectors: RITZ = Y*Z2;
   Tstart3 = MPI_Wtime()

   if (rank.eq.0) then
      do i = 1, nrhs
         print *, '>>>  ', W2(i)
      enddo
   endif

!!!!    C := alpha*op( A )*op( B ) + beta*C,
!!!!    RITZ = 1 * Y * Z2 + 0 * RITZ
!!!!    DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

    trans1='N';  trans2='N'
    alpha = 1.0d0;   beta = 0.0d0
    call dgemm(trans1,trans2,nloc,nrhs,nrhs,alpha,Y,nloc,Z2,ldz2,beta,RITZ,nloc)

!   do i = 1, nloc
!      do j = 1, nrhs
!         RITZ(i,j) = 0.0d0
!         do k = 1, nrhs
!            RITZ(i,j)=RITZ(i,j)+Y(i,k)*Z2(k,j)
!         enddo
!      enddo
!   enddo
   T_ritz_4 = T_ritz_4 + MPI_Wtime() - Tstart3

!....... Y = RITZ;
   do j = 1, nrhs
      do i = 1, nloc;   Y(i,j) = RITZ(i,j);   enddo
   enddo

   T_ritz = T_ritz + MPI_Wtime() - Tstart2
! -------------------------      
9090 continue   
   !-------------------
   ! move converged eigenvectors
!!   goto 9092

   Tstart2 = MPI_Wtime()

   call move_converged_eigenvectors(neqns,nzmax,ia,ja,a,Y,Yconv,AY,b,nloc,nrhs, &
            bandwidth,tracemin_tol, abs_res_tol, lam_tol,  &
            spikeparm,num_procs,rank,nconv,ndesired, &
            lam, resid_norm, Jvector, T_conv_1, T_conv_2, xtmp, xtmp_loc, min_res_norm,norma, &
            min_lam_curr, min_lam_prev, min_lam_rel_change)

   T_conv = T_conv + MPI_Wtime() - Tstart2

   if (debug .and. (rank .eq. 0)) print *, '      Min residual norm = ', min_res_norm
   if (nconv .ge. ndesired) goto 9080
9092 continue
!!-------------------
!! deflate
!!   goto 9093
!!   if (nconv .ge. 1) then
!!      Tstart2 = MPI_Wtime()
!!      call apply_deflation(nloc,nrhs,nconv,Y,Yconv,B_loc,B_glob, Bmat, ldBmat)
!!      T_defl = T_defl + MPI_Wtime() - Tstart2
!!   endif  !  if (nconv .ge. 1) then (deflate)
!!   T_other = T_other + MPI_Wtime() - Tstart
!!
!!9093 continue
   !-------------------
enddo !  do trmin_iters = 1, max_trmin_iters

9080 continue
!
!==================================================================================

! write the fiedler vector
allocate(rcounts(num_procs)) 
allocate(displs(num_procs)) 
do i =1,num_procs-1
  rcounts(i) = part1
  displs(i) = (i-1)*part1
end do
rcounts(num_procs) = part1_LAST
displs(num_procs) = (num_procs-1)*part1

!allocate(fiedler_vec(neqns)) 
call MPI_GATHERV(Yconv(1:nloc,2) , nloc ,MPI_DOUBLE_PRECISION,fiedler_vec, rcounts, displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,error)
!if ( rank .eq. 0 ) then 
!call KB07AD (fiedler_vec, neqns, perm)
!do i=1,neqns  
!write(10,*) fiedler_vec(i), perm(i) 
!end do 
!end if 

call form_resid_norms(neqns, nzmax, ia, ja, a, Yconv, AY, b, nloc, nconv, bandwidth, tol, &
                 spikeparm, num_procs, rank, lam, resid_norm, xtmp, xtmp_loc) 

do j = 1, nconv;   if (debug .and. (rank.eq.0)) print *, ' converged eigs : ', j, lam(j), resid_norm(j)/norma;   enddo


if (rank .eq. 0) print *, ' lam(2) / norma = ', lam(2)/norma

call form_resid_norms(neqns, nzmax, ia, ja, a, Y, AY, b, nloc, nrhs, bandwidth, tol, &
                 spikeparm, num_procs, rank, lam, resid_norm, xtmp, xtmp_loc) 

do j = 1, nrhs;   if (debug .and. (rank.eq.0)) print *, ' approx. eigs   : ', j, lam(j), resid_norm(j)/norma;   enddo


!write(40+rank,*) Yconv(:,2) 



!
!==================================================================================
!==================================================================================
!==================================================================================
!
   job = 3;
   if (rank .eq. 0) nzmax = ia(neqns+1)-1
   call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,error)
   !call SPIKEPARDISO(job, neqns, nzmax, ia, ja, a, b,nloc, nrhs,  bandwidth, tol,  spikeparm,  num_procs, rank );
   call MATVEC(job, neqns, nzmax, ia, ja, a, b,nloc, nrhs, num_procs, rank )
!
!==================================================================================
!==================================================================================
!==================================================================================
!
if (debug) then 
if (rank .eq. 0) then
   print *, ' '
   print *, ' =========================================== '
   print *, ' '
   print *, ' Performed ', trmin_iters, ' iterations of Trace Mininmization'
   print *, ' tracemin_tol = ', tracemin_tol
   print *, ' pcg_tol = ', pspike_tol
   print *, ' Number of desired eigenvalues = ', ndesired
   print *, ' nx = ', nx, ' N = ', neqns, ',   bandwidth = ', bandwidth

603 format(' spikeparams(1 to 8) = ', 8i3)

   print *, ' Number of MPI processes = ',num_procs
   if (inner_solver_flag .eq. 1) then
      print *, ' Inner solver = SPIKE-PARDSIO '
      write(6,603) spikeparm(1),spikeparm(2),spikeparm(3),spikeparm(4), &
                spikeparm(5),spikeparm(6),spikeparm(7),spikeparm(8) 
   elseif (inner_solver_flag .eq. 2) then
      print *, ' Inner solver = CG '
   endif

   print *, ' '
   print *, ' T_preprocess1 = ', T_preprocess1
   print *, ' T_preprocess2 = ', T_preprocess2
   print *, ' T_solve       = ', T_solve
   print *, '    T_cgmatvec = ', T_cgmatvec
   print *, '    T_cg_allred= ', T_cg_allred
   print *, ' T_other       = ', T_other

   print *, '      T_trmin  = ', T_trmin
   print *, '          [1]  = ', T_trmin_1
   print *, '          [2]  = ', T_trmin_2
   print *, '          [3]  = ', T_trmin_3
   print *, '          [4]  = ', T_trmin_4
   print *, '      T_ritz   = ', T_ritz
   print *, '          [1]  = ', T_ritz_1
   print *, '          [2]  = ', T_ritz_2
   print *, '          [3]  = ', T_ritz_3
   print *, '          [4]  = ', T_ritz_4
   print *, '      T_conv   = ', T_conv
   print *, '          [1]  = ', T_conv_1
   print *, '          [2]  = ', T_conv_2
   print *, '      T_orthog = ', T_orthog
   print *, '      T_defl   = ', T_defl
   print *, '      sum      = ', T_trmin+T_ritz+T_conv+T_orthog+T_defl
   print *, ' Total iteration time = ', T_solve + T_other 

   print *, ' nconv         = ', nconv 
   print *, ' nrhs          = ', nrhs

endif ! if (rank .eq. 0) then
end if 

   job =999 ! clean up
   call MATVEC(job, neqns, nzmax, ia, ja, a, b,nloc, nrhs, num_procs, rank )


   deallocate(b)
   deallocate(Y)
   deallocate(W)
   deallocate(AY)
   deallocate(RITZ);    
   deallocate(WORK)
   deallocate(Yconv)
   deallocate(xtmp, xtmp_loc, lam, resid_norm, Jvector)
   deallocate(Schur, Schur_loc_buf, Schur_glob_buf)
   deallocate(lu_WORK)
   deallocate(lu_IPIV)
   deallocate(W2, Z2, work2)
   deallocate(H_glob, H_loc)
   deallocate(B_loc, B_glob, Bmat)
   deallocate(rcounts)
   deallocate(displs) 


 
end subroutine FIEDLER 

!=======================================================
      subroutine get_nnzA_laplacian3d(nx, N, nnzA)
      implicit none

      integer nx, ny, nz, nnzA, N

      integer i, j, k, index, row

      ny=nx
      nz=nx
      N = nx*ny*nz

      nnzA=0
      
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (k .ne. 1) nnzA = nnzA + 1
               if (j .ne. 1) nnzA = nnzA + 1
               if (i .ne. 1) nnzA = nnzA + 1
               nnzA = nnzA + 1
            enddo
         enddo
      enddo
         
      return
    end subroutine get_nnzA_laplacian3d


subroutine form_laplacian3d(nx, N, nnzA, a, ia, ja)
  implicit none
  
  double precision A(1), valdiag, val
  integer IA(1), JA(1), nx, ny, nz, nnzA, N
  
  integer i, j, k, index, row
  character*50 str
  
  ny=nx;      nz=nx;      valdiag = 6.0d0;      val = -1.0d0
  
  row = 1;      index = 1
  ia(row) = 1  
  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           a(index) = valdiag;     ja(index) = row 
           index = index + 1

           if (i .ne. nx) then
              a(index) = val;      ja(index) = row + 1
              index = index + 1
           endif

           if (j .ne. ny) then
              a(index) = val;      ja(index) = row + nx  
              index = index + 1
           endif

           if (k .ne. nz) then
              a(index) = val;      ja(index) = row + nx*ny  
              index = index + 1
           endif
           
           row = row + 1
           ia(row) = index
        enddo
     enddo
  enddo
  
  return
end subroutine form_laplacian3d
!=======================================================
double precision function myrand(seed)
implicit none
double precision seed

integer i, iseed

iseed = nint(seed*1024)
i = mod(577*iseed+1,1024)
seed = dfloat(i)/1024.

myrand = seed

return
end function myrand
!=======================================================

!=======================================================


!=======================================================


!=======================================================
subroutine par_orthogonal_basis(nloc, nrhs, Q, tmp, tmp_loc)
  use mpi
  implicit none
  integer nloc, nrhs
  double precision Q(nloc,1)
  integer i, j, k, ierr
  double precision tmp(1), tmp_loc(1), ddot, dnrm2, tmpnrm, dsqrt

!!!  call par_dnrm2_single(nloc, Q(1,1), tmpnrm)  
  tmp_loc(1) = dnrm2(nloc,Q(1,1),1)**2

  call MPI_Allreduce(tmp_loc(1), tmp(1), 1, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)
  tmp(1) = dsqrt(tmp(1))

  do k = 1, nloc;   Q(k,1)=Q(k,1)/tmp(1);   enddo

  do i = 2, nrhs
     do j = 1, i-1
         tmp_loc(j) = ddot(nloc, Q(1,i), 1, Q(1,j), 1)
     enddo
     ! note that this is classical GS, to have fewer allreduces
     call MPI_Allreduce(tmp_loc, tmp, i-1, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)
     do j = 1, i-1
        do k = 1, nloc
           Q(k,i) = Q(k,i) - tmp(j)*Q(k,j)
        enddo
     enddo
!!!     call par_dnrm2_single(nloc, Q(1,i), tmpnrm)  
     tmp_loc(1) = dnrm2(nloc,Q(1,i),1)**2

     call MPI_Allreduce(tmp_loc(1), tmp(1), 1, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)
     tmp(1) = dsqrt(tmp(1))

     do k = 1, nloc;   Q(k,i) = Q(k,i)/tmp(1);   enddo
  enddo
 
  return
end subroutine par_orthogonal_basis
!=======================================================
subroutine par_dnrm2_single(neqns_loc, Y, xnrm)
  use mpi 
  implicit none
  integer neqns_loc
  double precision Y(1), xnrm

  integer j, ierr
  double precision dnrm2, xnrm_loc

  xnrm_loc = dnrm2(neqns_loc,Y(1),1)

  call MPI_Allreduce(xnrm_loc, xnrm, 1, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)

  return
end subroutine par_dnrm2_single
!==========================================================         
subroutine par_dnrm2_multiple(neqns_loc, nrhs, Y, xnrm, xnrm_loc)
  use mpi 

  implicit none

  integer neqns_loc, nrhs
  double precision Y(neqns_loc,1), xnrm(1), xnrm_loc(1)

  integer j, ierr
  double precision dnrm2

  do j = 1, nrhs
     xnrm_loc(j) = dnrm2(neqns_loc,Y(1,j),1)
  enddo

  call MPI_Allreduce(xnrm_loc, xnrm, nrhs , MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)

  return
end subroutine par_dnrm2_multiple
!=======================================================
subroutine par_ddot_single(neqns_loc, X, Y, xddot)
  use mpi 
  implicit none
  integer neqns_loc
  double precision X(1), Y(1), xddot

  integer j, ierr
  double precision ddot, xddot_loc

  xddot_loc = ddot(neqns_loc,X(1),1,Y(1),1)

  call MPI_Allreduce(xddot_loc, xddot, 1, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)

  return
end subroutine par_ddot_single
!=======================================================
subroutine form_resid_norms(neqns, nzmax, ia, ja, a, Y, AY, b, nloc, nrhs, bandwidth, tol, &
                 spikeparm, num_procs, rank, lam, resid_norm, xtmp, xtmp_loc) 
use mpi
implicit none 

integer neqns, nzmax, ia(1), ja(1), nloc, nrhs, bandwidth
integer  spikeparm(1), num_procs, rank
double precision a(1), Y(nloc,1), AY(nloc,1), b(nloc,1), tol, lam(1), resid_norm(1)
double precision  xtmp(1), xtmp_loc(1), ddot, dnrm2, dsqrt

integer i, j, job, ierr
integer flag

flag = 1

if (flag.eq.1) then
!   AY = Y
   do j = 1, nrhs
      do i = 1, nloc; AY(i,j) = Y(i,j); enddo
   enddo

   job = 55; !! JOB = 2 for solve 
   if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
   call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!   call SPIKEPARDISO(job, neqns, nzmax, ia, ja, a, AY, nloc, nrhs,  bandwidth, tol,  &
!        spikeparm, num_procs, rank );
   call MATVEC(job, neqns, nzmax, ia, ja, a, AY,nloc, nrhs, num_procs, rank )

   do j = 1, nrhs;
      xtmp_loc(j) = ddot(nloc, Y(1,j), 1, AY(1,j), 1)  !!!  lam(j)); 
   enddo

   call MPI_Allreduce(xtmp_loc, lam, nrhs, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)

   do j = 1, nrhs;
      do i = 1, nloc;   AY(i,j) = AY(i,j) - lam(j) * Y(i,j);   enddo

      !xtmp_loc(j) = dnrm2(nloc, AY(1,j), 1)**2
       xtmp_loc(j) = maxval(abs(AY(1:nloc,j))) 
   enddo

!   call MPI_Allreduce(xtmp_loc, resid_norm, nrhs, MPI_DOUBLE_PRECISION,  &
!                   MPI_SUM, MPI_COMM_WORLD, ierr)
!   do j = 1, nrhs
!      resid_norm(j) = dsqrt(resid_norm(j))
!   enddo
call MPI_Allreduce(xtmp_loc, resid_norm, nrhs, MPI_DOUBLE_PRECISION,  &
      MPI_MAX, MPI_COMM_WORLD, ierr)


else if (flag .eq. 2) then
   do j = 1, nrhs
      do i = 1, nloc; b(i,1) = Y(i,j); enddo
      job = 55; !! JOB = 2 for solve 
      if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
      call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!      call SPIKEPARDISO(job, neqns, nzmax, ia, ja, a, b, nloc, nrhs,  bandwidth, tol,  &
!                        spikeparm, num_procs, rank );
      call MATVEC(job, neqns, nzmax, ia, ja, a, b,nloc, nrhs, num_procs, rank )
      do i = 1, nloc; AY(i,j) = b(i,1); enddo
   enddo

!     re-write this to have 2 allreduces, not 2*nrhs
   do j = 1, nrhs
      call par_ddot_single(nloc, Y(1,j), AY(1,j), lam(j))
      do i = 1, nloc;   AY(i,j) = AY(i,j) - lam(j) * Y(i,j);   enddo
      call par_dnrm2_single(nloc, AY(1,j), resid_norm(j))
   enddo
else
   if (rank .eq. 0) print *, ' Wrong flag value in form_resid_norms '
endif

return
end subroutine form_resid_norms
!=======================================================
!    call move_converged_eigenvectors(neqns,nzmax,ia,ja,a,Y,Yconv,AY,b,nloc,nrhs, &
!                     bandwidth,tracemin_tol,spikeparm,num_procs,rank,nconv,ndesired, &
!                     lam, resid_norm, Jvector)

subroutine move_converged_eigenvectors(neqns,nzmax,ia,ja,a,Y,Yconv,AY,b,nloc,nrhs, &
                     bandwidth,tracemin_tol, abs_res_tol, lam_tol,  &
                     spikeparm,num_procs,rank,nconv,ndesired, &
                     lam, resid_norm, Jvector, T_conv_1, T_conv_2, xtmp, xtmp_loc, min_res_norm,norma, &
                     min_lam_curr, min_lam_prev, min_lam_rel_change)
  use mpi
  implicit none
  integer neqns, nzmax, nloc, nrhs, nconv, ndesired, ia(1), ja(1), trmin_iters2
  integer num_procs, rank
  double precision a(1), Y(nloc,1), Yconv(nloc,1), AY(nloc,1), b(nloc,1), lam(1), resid_norm(1)
  double precision tracemin_tol, xtmp(1), xtmp_loc(1), min_res_norm
  integer i, jj, j2, k, Jvector(nrhs), Jindex, spikeparm(1), bandwidth, ierr
  double precision rjjnorm, xjjnorm, dnrm2, ddot, dsqrt, T_conv_1, T_conv_2, Tstart
  double precision lam_min , norma 
  integer jmin
  double precision  min_lam_curr, min_lam_prev, min_lam_rel_change
  double precision  abs_res_tol, lam_tol

  Jindex = 0

  Tstart = MPI_Wtime()
  call form_resid_norms(neqns, nzmax, ia, ja, a, Y, AY, b, nloc, nrhs, bandwidth, tracemin_tol, &
       spikeparm, num_procs, rank, lam, resid_norm, xtmp, xtmp_loc) 
  T_conv_1 = T_conv_1 + MPI_Wtime() - Tstart

  min_res_norm = 1.0d10
  Tstart = MPI_Wtime() 

  lam_min = 1.d10
  jmin = 1
  do jj = 1, nrhs
     if (lam(jj) .lt. lam_min) then
        jmin = jj
        lam_min = lam(jj)
     endif
  enddo

  min_lam_curr = lam_min
  min_lam_rel_change = abs( (min_lam_curr - min_lam_prev) / min_lam_curr );
  min_lam_prev = min_lam_curr;
  if (rank.eq.0) print *, ' min_lam_curr       = ', min_lam_curr
  if (rank.eq.0) print *, ' min_lam_rel_change = ', min_lam_rel_change

  min_res_norm = resid_norm(jmin)/norma
  if (rank.eq.0) print *, ' min_res_norm       = ', min_res_norm
  if (rank.eq.0) print *, ' abs_resid          = ', resid_norm(jmin)

!!!  These are passed in now
!!!  abs_res_tol = 1e-5
!!!  lam_tol = 0.25

  do jj = 1, nrhs
!!!!!!  min_res_norm = min(min_res_norm, resid_norm(jj)/norma)     
     if (jj .eq. jmin .and. &
!
          ((resid_norm(jj)/norma  .le. tracemin_tol .and.  &
             resid_norm(jj) .le. abs_res_tol) &  
             .or. &
           (resid_norm(jj)/norma  .le. tracemin_tol .and.  &
             min_lam_rel_change .le. lam_tol) &
             .or. &
           (resid_norm(jj) .le. abs_res_tol .and.  &
             min_lam_rel_change .le. lam_tol)))   then
!
        nconv = nconv + 1
        !  is this necessary?
!!!        call par_dnrm2_single(nloc, Y(1,jj), xjjnorm)
           xtmp_loc(1) = dnrm2(nloc,Y(1,jj),1)**2

           call MPI_Allreduce(xtmp_loc(1), xtmp(1), 1, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)
           xjjnorm = dsqrt(xtmp(1))


        do i = 1, nloc;   Yconv(i,nconv) = Y(i,jj)/xjjnorm;   enddo
        if (rank .eq. 0) print *, ' ====== Eigenvector converged === ', lam(jj),  resid_norm(jj)/norma
        !  if (nconv .ge. ndesired) goto 999
     else
           Jindex = Jindex + 1;   Jvector(Jindex) = jj
     endif
  enddo  !!    do jj = 1, nrhs

  T_conv_2 = T_conv_2 + MPI_Wtime() - Tstart

!  Shift columns of X left to make them contiguous
999  continue
     
  if (Jindex .lt. nrhs) then
     do j2 = 1, Jindex;   
        do i = 1, nloc;   Y(i,j2) = Y(i,Jvector(j2));   enddo
     enddo
  endif

  nrhs = Jindex   !  Modify nrhs:
       
  return
end subroutine move_converged_eigenvectors
!=======================================================
subroutine apply_deflation(nloc,nrhs,nconv,Y,Yconv,B_loc,B_glob, Bmat, ldBmat)
use mpi
implicit none

integer nloc, nrhs, nconv, ldBmat
double precision Y(nloc,1), Yconv(nloc,1), B_loc(1), B_glob(1), Bmat(ldBmat,1), ddot
character*1  trans1, trans2
double precision alpha, beta
integer i, j, k, ierr

k = 1
do i = 1, nconv
   do j = 1, nrhs
!     B(i,j)=ddot(n,Yconv(1,i),1,Y(1,j),1)
      B_loc(k) = ddot(nloc,Yconv(1,i),1,Y(1,j),1)
      k = k+1
   enddo
enddo

call MPI_Allreduce(B_loc, B_glob, k-1, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)

k = 1
do i = 1, nconv
   do j = 1, nrhs
      Bmat(i,j) = B_glob(k)
      k = k+1
   enddo
enddo

!!!!    C := alpha*op( A )*op( B ) + beta*C,
!!!!    Y = Y - Yconv * Bmat
!!!!    DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

 trans1='N';  trans2='N'
 alpha = -1.0d0;   beta = 1.0d0
 call dgemm(trans1,trans2,nloc,nrhs,nconv,alpha,Yconv,nloc,Bmat,ldBmat,beta,Y,nloc)

!do i = 1, nloc
!   do j = 1, nrhs
!!!!!!!      Y(i,j) = Y(i,j) 
!      do k = 1, nconv
!         Y(i,j) = Y(i,j) - Yconv(i,k)*Bmat(k,j)
!      enddo
!   enddo
!enddo

return
end subroutine apply_deflation
!=======================================================
subroutine conjugate_gradient(neqns, nzmax, ia, ja, a, B, X, nloc, nrhs, bandwidth, tol, &
                 spikeparm, num_procs, rank, xtmp, xtmp_loc, &
                 P, Y, W, R, iter, maxit, T_cgmatvec, T_cg_allred) 
  use mpi
  implicit none 

  integer neqns, nzmax, ia(1), ja(1), nloc, nrhs, bandwidth, iter, maxit
  integer  spikeparm(1), num_procs, rank
  double precision a(1), tol, B(nloc,1), X(nloc,1) 
  double precision P(nloc,1), Y(nloc,1), W(nloc,1), R(nloc,1)
  double precision  xtmp(1), xtmp_loc(1), ddot, dnrm2, dsqrt, T_cgmatvec, T_cg_allred, Tstart 
  
  integer i, j, job, ierr

  double precision, dimension(:), allocatable :: prn, prn0, alpha, beta, mu, rho, rho_old
!=============================================================
  allocate(prn(nrhs), prn0(nrhs), alpha(nrhs), beta(nrhs), mu(nrhs), rho(nrhs), rho_old(nrhs))
!=============================================================
!!!X = zeros(size(B));
!!![NN,k]=size(B);
!!!R = B;   %Y = A*X; R = B - Y;
  do j = 1, nrhs
     do i = 1, nloc;  X(i,j) = 0.0d0; R(i,j) = B(i,j);  Y(i,j) = R(i,j); enddo
  enddo
!!!
!!!Y=A*R - shift*R; 
!  do j = 1, nrhs
!     do i = 1, nloc; Y(i,j) = R(i,j); enddo
!  enddo

  Tstart = MPI_Wtime()
  job = 55; !! JOB = 2 for solve 
  if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
  call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!  call SPIKEPARDISO(job, neqns, nzmax, ia, ja, a, Y, nloc, nrhs,  bandwidth, tol,  &
!       spikeparm, num_procs, rank );
  call MATVEC(job, neqns, nzmax, ia, ja, a, Y,nloc, nrhs, num_procs, rank )
  T_cgmatvec = T_cgmatvec + MPI_Wtime() - Tstart
!!!--------------------
!!!beta = zeros(k,1); alpha = beta; prn = beta; prn0 = beta;
!!!rho = beta; mu = beta; rho_old=ones(k,1); 

  do j = 1, nrhs
     alpha(j) = 0.0d0;   beta(j) = 0.0d0;   prn(j) = 0.0d0;   prn0(j) = 0.0d0;   
     rho(j) = 0.0d0;   rho_old(j) = 1.0d0;   mu(j) = 0.0d0;   
  enddo

!!!--------------------
!!!for i = 1:k, prn(i) = norm(R(:,i));        end % =====
!!!prn0 = prn;

 do j = 1, nrhs;
    xtmp_loc(j) = dnrm2(nloc, R(1,j), 1)**2
 enddo

 Tstart = MPI_Wtime()
 call MPI_Allreduce(xtmp_loc, prn, nrhs, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)
 T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart

 do j = 1, nrhs
    prn(j) = dsqrt(prn(j));  prn0(j) = prn(j)
  enddo
!!!--------------------
!!!W=zeros(size(X)); P=W; 
  do j = 1, nrhs
     do i = 1, nloc; W(i,j) = 0.0d0;  P(i,j) = 0.0d0; enddo
  enddo
!!!--------------------
!!!iter = 0;
!!!while (iter<maxit && norm(prn) > rtol* norm(prn0))
  do iter = 1, maxit
     !!!--------------------
     !!!   for i = 1:k, rho(i) = R(:,i)'*Y(:,i);   end % =====
     do j = 1, nrhs
        xtmp_loc(j) = ddot(nloc,R(1,j),1,Y(1,j),1)
     enddo

     Tstart = MPI_Wtime()
     call MPI_Allreduce(xtmp_loc, rho, nrhs, MPI_DOUBLE_PRECISION,  &
                        MPI_SUM, MPI_COMM_WORLD, ierr)
     T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart
     !!!--------------------
     !!!   for i = 1:k, beta(i) = rho(i) / rho_old(i); 
     !!!                rho_old(i) = rho(i);       end % =====
     do j = 1, nrhs
        beta(j) = rho(j) / rho_old(j); rho_old(j) = rho(j)
     enddo
     !!!--------------------
     !!!   P = R + P * diag(beta);
     !!!   W = Y + W * diag(beta);
!$OMP PARALLEL DO PRIVATE(i,j)
     do j = 1, nrhs
        do i = 1, nloc
           P(i,j) = R(i,j) + beta(j) * P(i,j); 
           W(i,j) = Y(i,j) + beta(j) * W(i,j); 
        enddo
     enddo
!$OMP PARALLEL DO PRIVATE(i,j)
     !!!--------------------
     !!!   for i = 1:k, mu(i) = W(:,i)'*W(:,i);    end % =====
     do j = 1, nrhs
        xtmp_loc(j) = ddot(nloc,W(1,j),1,W(1,j),1)
     enddo

     Tstart = MPI_Wtime()
     call MPI_Allreduce(xtmp_loc, mu, nrhs, MPI_DOUBLE_PRECISION,  &
                        MPI_SUM, MPI_COMM_WORLD, ierr)
     T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart
     !!!--------------------
     !!!   for i = 1:k, alpha(i) = rho(i)/mu(i);   end % =====u;
     do j = 1, nrhs; alpha(j) = rho(j)/mu(j); enddo
     !!!--------------------
     !!!   X = X + P * diag(alpha);
     !!!   R = R - W * diag(alpha);
!$OMP PARALLEL DO PRIVATE(i,j)
     do j = 1, nrhs
        do i = 1, nloc
           X(i,j) = X(i,j) + alpha(j) * P(i,j)
           R(i,j) = R(i,j) - alpha(j) * W(i,j)
           Y(i,j) = R(i,j)
         enddo
     enddo
!$OMP END PARALLEL DO
     !!!--------------------
     !!!   Y = A*R - shift*R; 
!     do j = 1, nrhs
!        do i = 1, nloc; Y(i,j) = R(i,j); enddo
!     enddo

     Tstart = MPI_Wtime()
     job = 55; !! JOB = 2 for solve 
     if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
     call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!     call SPIKEPARDISO(job, neqns, nzmax, ia, ja, a, Y, nloc, nrhs,  bandwidth, tol,  &
!            spikeparm, num_procs, rank );
     call MATVEC(job, neqns, nzmax, ia, ja, a, Y,nloc, nrhs, num_procs, rank )
     T_cgmatvec = T_cgmatvec + MPI_Wtime() - Tstart
     !!!--------------------
     !!!   for i = 1:k, prn(i) = norm(R(:,i));     end % =====
!     do j = 1, nrhs;
!        xtmp_loc(j) = dnrm2(nloc, R(1,j), 1)**2
!     enddo
!
!     Tstart = MPI_Wtime()
!     call MPI_Allreduce(xtmp_loc, prn, nrhs, MPI_DOUBLE_PRECISION,  &
!                   MPI_SUM, MPI_COMM_WORLD, ierr)
!     T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart
!     do j = 1, nrhs
!        prn(j) = dsqrt(prn(j));
!     enddo
     !!!--------------------
     !!!   iter = iter+1;
  !!!end
  enddo

return
end subroutine conjugate_gradient






!=======================================================
subroutine conjugate_gradient_mod(neqns, nzmax, ia, ja, a, B, X, nloc, nrhs, bandwidth, tol, &
                 spikeparm, num_procs, rank, xtmp, xtmp_loc, &
                 P, Y, W, R, iter, maxit, T_cgmatvec, T_cg_allred) 
  use mpi
  implicit none 

  integer neqns, nzmax, ia(1), ja(1), nloc, nrhs, bandwidth, iter, maxit
  integer  spikeparm(1), num_procs, rank
  double precision a(1), tol, B(nloc,1), X(nloc,1) 
  double precision P(nloc,1), Y(nloc,1), W(nloc,1), R(nloc,1)
  double precision  xtmp(1), xtmp_loc(1), ddot, dnrm2, dsqrt, T_cgmatvec, T_cg_allred, Tstart 
  
  integer i, j, job, ierr

  double precision, dimension(:), allocatable :: prn, prn0, alpha, beta, mu, rho, rho_old
  double precision, dimension(:), allocatable :: gamma, delta, x3_loc, x3_glob
!=============================================================
  allocate(prn(nrhs), prn0(nrhs), alpha(nrhs), beta(nrhs), mu(nrhs), rho(nrhs), rho_old(nrhs))
  allocate(gamma(nrhs), delta(nrhs),  x3_loc(3*nrhs), x3_glob(3*nrhs))
!=============================================================
!!!X = zeros(size(B));
!!![NN,k]=size(B);
!!!R = B;   %Y = A*X; R = B - Y;
  do j = 1, nrhs
     do i = 1, nloc;  X(i,j) = 0.0d0; R(i,j) = B(i,j);  Y(i,j) = R(i,j); enddo
  enddo
!!!
!!!Y=A*R - shift*R; 
!  do j = 1, nrhs
!     do i = 1, nloc; Y(i,j) = R(i,j); enddo
!  enddo

  Tstart = MPI_Wtime()
  job = 55; !! JOB = 2 for solve 
  if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
  call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!  write(10+rank,*) Y 
  call MATVEC(job, neqns, nzmax, ia, ja, a, Y,nloc, nrhs, num_procs, rank )
!  write(20+rank,*) Y 
!   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!   stop 


  T_cgmatvec = T_cgmatvec + MPI_Wtime() - Tstart
!!!--------------------
!!!beta = zeros(k,1); alpha = beta; prn = beta; prn0 = beta;
!!!rho = beta; mu = beta; rho_old=ones(k,1); 

  do j = 1, nrhs
     alpha(j) = 0.0d0;   beta(j) = 0.0d0;   prn(j) = 0.0d0;   prn0(j) = 0.0d0;   
     rho(j) = 0.0d0;   rho_old(j) = 1.0d0;   mu(j) = 0.0d0;   
  enddo

!!!--------------------
!!!for i = 1:k, prn(i) = norm(R(:,i));        end % =====
!!!prn0 = prn;

 do j = 1, nrhs;
    xtmp_loc(j) = dnrm2(nloc, R(1,j), 1)**2
 enddo

 Tstart = MPI_Wtime()
 call MPI_Allreduce(xtmp_loc, prn, nrhs, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)
 T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart

 do j = 1, nrhs
    prn(j) = dsqrt(prn(j));  prn0(j) = prn(j)
  enddo
!!!--------------------
!!!W=zeros(size(X)); P=W; 
  do j = 1, nrhs
     do i = 1, nloc; W(i,j) = 0.0d0;  P(i,j) = 0.0d0; enddo
  enddo
!!!--------------------
!!!iter = 0;
!!!while (iter<maxit && norm(prn) > rtol* norm(prn0))
  do iter = 1, maxit
     !!!--------------------
     !!!   for i = 1:k, rho(i) = R(:,i)'*Y(:,i);   end % =====
     do j = 1, nrhs
        x3_loc(j       ) = ddot(nloc,R(1,j),1,Y(1,j),1) ! rho
        x3_loc(j+  nrhs) = ddot(nloc,Y(1,j),1,Y(1,j),1) ! gamma
        x3_loc(j+2*nrhs) = ddot(nloc,W(1,j),1,Y(1,j),1) ! delta
     enddo

     Tstart = MPI_Wtime()
!     call MPI_Allreduce(x3_loc, x3_glob, 3*nrhs, MPI_DOUBLE_PRECISION,  &
!                        MPI_SUM, MPI_COMM_WORLD, ierr)


     call MPI_Reduce(x3_loc, x3_glob, 3*nrhs, MPI_DOUBLE_PRECISION,  &
                        MPI_SUM, 0, MPI_COMM_WORLD, ierr)

     call MPI_BCAST(x3_glob , 3*nrhs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
     

     T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart

     do j = 1, nrhs
        rho(j)    = x3_glob(j       ) !! ddot(nloc,R(1,j),1,Y(1,j),1) ! rho
        gamma(j)  = x3_glob(j+  nrhs) !! ddot(nloc,Y(1,j),1,Y(1,j),1) ! gamma
        delta(j)  = x3_glob(j+2*nrhs) !! ddot(nloc,W(1,j),1,Y(1,j),1) ! delta
     enddo


     !!!--------------------
     !!!   for i = 1:k, beta(i) = rho(i) / rho_old(i); 
     !!!                rho_old(i) = rho(i);       end % =====
     do j = 1, nrhs
        beta(j) = rho(j) / rho_old(j); rho_old(j) = rho(j)
     enddo
     !!!--------------------
     !!!   P = R + P * diag(beta);
     !!!   W = Y + W * diag(beta);
!$OMP PARALLEL DO PRIVATE(i,j)
     do j = 1, nrhs
        do i = 1, nloc
           P(i,j) = R(i,j) + beta(j) * P(i,j); 
           W(i,j) = Y(i,j) + beta(j) * W(i,j); 
        enddo
     enddo
!$OMP PARALLEL DO PRIVATE(i,j)
     !!!--------------------
     !!!   for i = 1:k, mu(i) = W(:,i)'*W(:,i);    end % =====
     do j = 1, nrhs
!!!!!!        xtmp_loc(j) = ddot(nloc,W(1,j),1,W(1,j),1)
     enddo

     Tstart = MPI_Wtime()
!!!!!!     call MPI_Allreduce(xtmp_loc, mu, nrhs, MPI_DOUBLE_PRECISION,  &
!!!!!!                        MPI_SUM, MPI_COMM_WORLD, ierr)
     T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart
     !!!--------------------
     !!!   for i = 1:k, alpha(i) = rho(i)/mu(i);   end % =====u;
     do j = 1, nrhs
        mu(j) = gamma(j) + beta(j) * (2.d0*delta(j) + beta(j)*mu(j)) 
        alpha(j) = rho(j)/mu(j); 
     enddo
     !!!--------------------
     !!!   X = X + P * diag(alpha);
     !!!   R = R - W * diag(alpha);
!$OMP PARALLEL DO PRIVATE(i,j)
     do j = 1, nrhs
        do i = 1, nloc
           X(i,j) = X(i,j) + alpha(j) * P(i,j)
           R(i,j) = R(i,j) - alpha(j) * W(i,j)
           Y(i,j) = R(i,j)
         enddo
     enddo
!$OMP END PARALLEL DO
     !!!--------------------
     !!!   Y = A*R - shift*R; 
!     do j = 1, nrhs
!        do i = 1, nloc; Y(i,j) = R(i,j); enddo
!     enddo

     Tstart = MPI_Wtime()
     job = 55; !! JOB = 2 for solve 
     if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
     call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!     call SPIKEPARDISO(job, neqns, nzmax, ia, ja, a, Y, nloc, nrhs,  bandwidth, tol,  &
!            spikeparm, num_procs, rank );
     call MATVEC(job, neqns, nzmax, ia, ja, a, Y,nloc, nrhs, num_procs, rank )
     T_cgmatvec = T_cgmatvec + MPI_Wtime() - Tstart
     !!!--------------------
     !!!   for i = 1:k, prn(i) = norm(R(:,i));     end % =====
!     do j = 1, nrhs;
!        xtmp_loc(j) = dnrm2(nloc, R(1,j), 1)**2
!     enddo
!
!     Tstart = MPI_Wtime()
!     call MPI_Allreduce(xtmp_loc, prn, nrhs, MPI_DOUBLE_PRECISION,  &
!                   MPI_SUM, MPI_COMM_WORLD, ierr)
!     T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart
!     do j = 1, nrhs
!        prn(j) = dsqrt(prn(j));
!     enddo
     !!!--------------------
     !!!   iter = iter+1;
  !!!end
  enddo

return
end subroutine conjugate_gradient_mod
!=======================================================


subroutine pcr(prec, neqns, nzmax, ia, ja, a, B, X, nloc, nrhs, bandwidth, tol, &
                 spikeparm, num_procs, rank, xtmp, xtmp_loc, &
                 P, Y, W, R, iter, maxit, T_cgmatvec, T_cg_allred)
  use mpi
  implicit none
  logical prec 
  integer neqns, nzmax, ia(1), ja(1), nloc, nrhs, bandwidth, iter, maxit
  integer  spikeparm(1), num_procs, rank
  double precision a(1), tol, B(nloc,1), X(nloc,1)
  double precision P(nloc,1), Y(nloc,1), W(nloc,1), R(nloc,1)
  double precision  xtmp(1), xtmp_loc(1), ddot, dnrm2, dsqrt, T_cgmatvec, T_cg_allred, Tstart
   
  integer i, j, job, ierr

  double precision, dimension(:), allocatable :: prn, prn0, alpha, beta, mu, rho, rho_old

! murat 
  double precision, dimension(:,:) , allocatable :: MW 
! murat 

!=============================================================
  allocate(prn(nrhs), prn0(nrhs), alpha(nrhs), beta(nrhs), mu(nrhs), rho(nrhs), rho_old(nrhs))

! murat 
  allocate(MW(nloc, nrhs)) 
! murat 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
prec = .false. !  preconditioner still don't work for pcr 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=============================================================
!!!X = zeros(size(B));
!!![NN,k]=size(B);
!!!R = B;   %Y = A*X; R = B - Y;
  do j = 1, nrhs
     do i = 1, nloc;  X(i,j) = 0.0d0;MW(i,j) =0.0d0;  R(i,j) = B(i,j); enddo
  enddo


if (prec) then
      job = 88; !! JOB = 2 for solve
      if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
      call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MATVEC(job, neqns, nzmax, ia, ja, a, R,nloc, nrhs, num_procs, rank )
end if


  do j = 1, nrhs
     do i = 1, nloc; Y(i,j) = R(i,j); enddo
  enddo


!!!
!!!Y=A*R - shift*R;
!  do j = 1, nrhs
!     do i = 1, nloc; Y(i,j) = R(i,j); enddo
!  enddo

  Tstart = MPI_Wtime()
  job = 55; !! JOB = 2 for solve
  if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
  call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!  call SPIKEPARDISO(job, neqns, nzmax, ia, ja, a, Y, nloc, nrhs,  bandwidth, tol,  &
!       spikeparm, num_procs, rank );
  call MATVEC(job, neqns, nzmax, ia, ja, a, Y,nloc, nrhs, num_procs, rank )
  T_cgmatvec = T_cgmatvec + MPI_Wtime() - Tstart
!!!--------------------
!!!beta = zeros(k,1); alpha = beta; prn = beta; prn0 = beta;
!!!rho = beta; mu = beta; rho_old=ones(k,1);

  do j = 1, nrhs
     alpha(j) = 0.0d0;   beta(j) = 0.0d0;   prn(j) = 0.0d0;   prn0(j) = 0.0d0;
     rho(j) = 0.0d0;   rho_old(j) = 1.0d0;   mu(j) = 0.0d0;
  enddo

!!!--------------------
!!!for i = 1:k, prn(i) = norm(R(:,i));        end % =====
!!!prn0 = prn;

! do j = 1, nrhs;
!    xtmp_loc(j) = dnrm2(nloc, R(1,j), 1)**2
! enddo
!
! Tstart = MPI_Wtime()
! call MPI_Allreduce(xtmp_loc, prn0, nrhs, MPI_DOUBLE_PRECISION,  &
!                   MPI_SUM, MPI_COMM_WORLD, ierr)
! T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart
do j = 1, nrhs;
     xtmp_loc(j) = maxval(abs(R(1:nloc,j)))
enddo
call MPI_Allreduce(xtmp_loc, prn0, nrhs, MPI_DOUBLE_PRECISION,  &
      MPI_MAX, MPI_COMM_WORLD, ierr)


!!!--------------------
!!!W=zeros(size(X)); P=W;
  do j = 1, nrhs
     do i = 1, nloc; W(i,j) = 0.0d0;  P(i,j) = 0.0d0; enddo
  enddo
!!!--------------------
!!!iter = 0;
!!!while (iter<maxit && norm(prn) > rtol* norm(prn0))
  do iter = 1, maxit
     !!!--------------------
     !!!   for i = 1:k, rho(i) = R(:,i)'*Y(:,i);   end % =====
     do j = 1, nrhs
        xtmp_loc(j) = ddot(nloc,R(1,j),1,Y(1,j),1)
     enddo

     Tstart = MPI_Wtime()
     call MPI_Allreduce(xtmp_loc, rho, nrhs, MPI_DOUBLE_PRECISION,  &
                        MPI_SUM, MPI_COMM_WORLD, ierr)
     T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart
     !!!--------------------
     !!!   for i = 1:k, beta(i) = rho(i) / rho_old(i);
     !!!                rho_old(i) = rho(i);       end % =====
     do j = 1, nrhs
        beta(j) = rho(j) / rho_old(j); rho_old(j) = rho(j)
     enddo
     !!!--------------------
     !!!   P = R + P * diag(beta);
     !!!   W = Y + W * diag(beta);
!$OMP PARALLEL DO PRIVATE(i,j)
     do j = 1, nrhs
        do i = 1, nloc
           P(i,j) = R(i,j) + beta(j) * P(i,j);
           W(i,j) = Y(i,j) + beta(j) * W(i,j);
           MW(i,j) = W(i,j) ; 
        enddo
     enddo
!$OMP END PARALLEL DO 

  
if (prec) then 
      job = 88; !! JOB = 2 for solve
      if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
      call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MATVEC(job, neqns, nzmax, ia, ja, a, MW,nloc, nrhs, num_procs, rank )
end if 


     !!!--------------------
     !!!   for i = 1:k, mu(i) = W(:,i)'*W(:,i);    end % =====
     do j = 1, nrhs
        xtmp_loc(j) = ddot(nloc,W(1,j),1,MW(1,j),1)
     enddo

     Tstart = MPI_Wtime()
     call MPI_Allreduce(xtmp_loc, mu, nrhs, MPI_DOUBLE_PRECISION,  &
                        MPI_SUM, MPI_COMM_WORLD, ierr)
     T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart
     !!!--------------------
     !!!   for i = 1:k, alpha(i) = rho(i)/mu(i);   end % =====u;
     do j = 1, nrhs; alpha(j) = rho(j)/mu(j); enddo
     !!!--------------------
     !!!   X = X + P * diag(alpha);
     !!!   R = R - W * diag(alpha);
!$OMP PARALLEL DO PRIVATE(i,j)
     do j = 1, nrhs
        do i = 1, nloc
           X(i,j) = X(i,j) + alpha(j) * P(i,j)
           R(i,j) = R(i,j) - alpha(j) * MW(i,j)
           Y(i,j) = R(i,j)
         enddo
     enddo
!$OMP END PARALLEL DO
     !!!--------------------
     !!!   Y = A*R - shift*R;
!     do j = 1, nrhs
!        do i = 1, nloc; Y(i,j) = R(i,j); enddo
!     enddo

     Tstart = MPI_Wtime()
     job = 55; !! JOB = 2 for solve
     if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
     call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!     call SPIKEPARDISO(job, neqns, nzmax, ia, ja, a, Y, nloc, nrhs,  bandwidth, tol,  &
!            spikeparm, num_procs, rank );
     call MATVEC(job, neqns, nzmax, ia, ja, a, Y,nloc, nrhs, num_procs, rank )
     T_cgmatvec = T_cgmatvec + MPI_Wtime() - Tstart
     !!!--------------------
     !!!   for i = 1:k, prn(i) = norm(R(:,i));     end % =====
!     do j = 1, nrhs;
!        xtmp_loc(j) = dnrm2(nloc, R(1,j), 1)**2
!     enddo
!
!     Tstart = MPI_Wtime()
!     call MPI_Allreduce(xtmp_loc, prn, nrhs, MPI_DOUBLE_PRECISION,  &
!                   MPI_SUM, MPI_COMM_WORLD, ierr)
!     T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart
!     do j = 1, nrhs
!        prn(j) = dsqrt(prn(j));
!     enddo
     !!!--------------------
     !!!   iter = iter+1;
  !!!end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!Check the relative residual   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!do j = 1, nrhs;
!     xtmp_loc(j) = dnrm2(nloc, W(1,j), 1)**2
!enddo
!call MPI_Allreduce(xtmp_loc, prn, nrhs, MPI_DOUBLE_PRECISION,  &
!      MPI_SUM, MPI_COMM_WORLD, ierr)
!do j =1,nrhs
!   prn(j) = dsqrt(prn(j)) /dsqrt(prn0(j))
!end do

do j = 1, nrhs;
     xtmp_loc(j) =  maxval(abs(W(1:nloc,j)))
enddo
call MPI_Allreduce(xtmp_loc, prn, nrhs, MPI_DOUBLE_PRECISION,  &
      MPI_MAX, MPI_COMM_WORLD, ierr)
do j =1,nrhs
     prn(j) =  prn(j) / prn0(j)
end do


if (maxval(prn) .le. tol ) exit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  enddo

!if ((rank .eq. 0)) write(*,*) 'PCR (max) relative residual:',iter,  maxval(prn)

deallocate(MW) 
return
end subroutine pcr  


subroutine pcgsingular(prec, neqns, nzmax, ia, ja, a, B, X, nloc, nrhs, bandwidth, tol, &
                 spikeparm, num_procs, rank, xtmp, xtmp_loc, &
                 P, Y, W, R, iter, maxit, T_cgmatvec, T_cg_allred)
  use mpi
  implicit none
  logical prec
  integer neqns, nzmax, ia(1), ja(1), nloc, nrhs, bandwidth, iter, maxit
  integer  spikeparm(1), num_procs, rank
  double precision a(1), tol, B(nloc,1), X(nloc,1)
  double precision P(nloc,1), Y(nloc,1), W(nloc,1), R(nloc,1)
  double precision  xtmp(1), xtmp_loc(1), ddot, dnrm2, dsqrt, T_cgmatvec, T_cg_allred, Tstart

  integer i, j, job, ierr

  double precision, dimension(:), allocatable :: prn, prn0, alpha, beta, mu, rho, rho_old

!=============================================================
  allocate(prn(nrhs), prn0(nrhs), alpha(nrhs), beta(nrhs), mu(nrhs), rho(nrhs), rho_old(nrhs))


  do j = 1, nrhs
     do i = 1, nloc;  X(i,j) = 0.0d0; R(i,j) = B(i,j); enddo
  enddo
  do j = 1, nrhs
     do i = 1, nloc;  W(i,j) = R(i,j);   enddo
  enddo
  do j = 1, nrhs
     alpha(j) = 0.0d0;   beta(j) = 0.0d0;
  enddo

if (prec) then
  ! Z =   M \ R;
  job = 77; !! JOB = 2 for solve
  if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
  call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MATVEC(job, neqns, nzmax, ia, ja, a, W,nloc, nrhs, num_procs, rank )
end if

!do j = 1, nrhs;
!     xtmp_loc(j) = dnrm2(nloc, W(1,j), 1)**2
!enddo
!call MPI_Allreduce(xtmp_loc, prn0, nrhs, MPI_DOUBLE_PRECISION,  &
!      MPI_SUM, MPI_COMM_WORLD, ierr)
do j = 1, nrhs;
     xtmp_loc(j) = maxval(abs(W(1:nloc,j)))
enddo
call MPI_Allreduce(xtmp_loc, prn0, nrhs, MPI_DOUBLE_PRECISION,  &
      MPI_MAX, MPI_COMM_WORLD, ierr)

do iter = 1,maxit
  ! rho = R' * Z ;
  do j = 1, nrhs
     xtmp_loc(j) = ddot(nloc,R(1,j),1,W(1,j),1)
  enddo
  Tstart = MPI_Wtime()
  call MPI_Allreduce(xtmp_loc, rho, nrhs, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)
  T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart

if (iter .eq. 1) then
  ! P = Z ;
  do j = 1, nrhs
   do i = 1, nloc;  P(i,j) = W(i,j);   enddo
  enddo
else
   !! beta = rho / rho_old ;
   do j = 1, nrhs
     beta(j)  = rho(j)/rho_old(j);
   end do
   ! p = z + beta*p ;
!$OMP PARALLEL DO PRIVATE(i,j)
     do j = 1, nrhs
        do i = 1, nloc
            P(i,j) = W(i,j) + beta(j) * P(i,j);
        enddo
     enddo
!$OMP END PARALLEL DO
end if
do j=1,nrhs
   rho_old(j) = rho(j)
end do

  ! Ap = A*p;
!$OMP PARALLEL DO PRIVATE(i,j)
     do j = 1, nrhs
        do i = 1, nloc
            Y(i,j) = P(i,j);
        enddo
     enddo
!$OMP END PARALLEL DO
  Tstart = MPI_Wtime()
  job = 44; !! JOB = 2 for solve
  if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
  call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MATVEC(job, neqns, nzmax, ia, ja, a, Y,nloc, nrhs, num_procs, rank )
  T_cgmatvec = T_cgmatvec + MPI_Wtime() - Tstart

  ! mu = P' * AP
  do j = 1, nrhs
     xtmp_loc(j) = ddot(nloc,P(1,j),1,Y(1,j),1)
  enddo
  Tstart = MPI_Wtime()
  call MPI_Allreduce(xtmp_loc, mu, nrhs, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)
  T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart
  ! alpha = rho / mu
  ! alpha = rho / mu
  do j = 1, nrhs
     alpha(j)  = rho(j)/mu(j);
  end do

!$OMP PARALLEL DO PRIVATE(i,j)
     do j = 1, nrhs
        do i = 1, nloc
           X(i,j) = X(i,j) + alpha(j) * P(i,j)
           R(i,j) = R(i,j) - alpha(j) * Y(i,j)
           W(i,j) = R(i,j)
         enddo
     enddo
!$OMP END PARALLEL DO


! Z = D\R;
if (prec) then
  job = 77; !! JOB = 2 for solve
  if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
  call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MATVEC(job, neqns, nzmax, ia, ja, a, W,nloc, nrhs, num_procs, rank )
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!Check the relative residual   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do j = 1, nrhs;
     xtmp_loc(j) =  maxval(abs(W(1:nloc,j)))
enddo
call MPI_Allreduce(xtmp_loc, prn, nrhs, MPI_DOUBLE_PRECISION,  &
      MPI_MAX, MPI_COMM_WORLD, ierr)
do j =1,nrhs
     prn(j) =  prn(j) / prn0(j)
end do

!do j = 1, nrhs;
!     xtmp_loc(j) = dnrm2(nloc, W(1,j), 1)**2
!enddo
!call MPI_Allreduce(xtmp_loc, prn, nrhs, MPI_DOUBLE_PRECISION,  &
!      MPI_SUM, MPI_COMM_WORLD, ierr)
!do j =1,nrhs
!   prn(j) = sqrt(prn(j)) /sqrt(prn0(j))
!end do
if (maxval(prn) .le. tol ) exit



!! if not converging in 10 iteration exit 
!if ((iter .eq. 10) .and. (maxval(prn) .ge. 1.0 )) exit 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end do

if (rank .eq. 0 ) write(*,*) 'CG (max) relative residual:', iter, maxval(prn)
!if (debug .and. (rank .eq. 0 )) write(*,*) 'CG (max) relative residual:', iter, maxval(prn)
return
end subroutine pcgsingular


subroutine pcg( prec, neqns, nzmax, ia, ja, a, B, X, nloc, nrhs, bandwidth, tol, &
                 spikeparm, num_procs, rank, xtmp, xtmp_loc, &
                 P, Y, W, R, iter, maxit, T_cgmatvec, T_cg_allred)
  use mpi
  implicit none
  logical prec
  integer neqns, nzmax, ia(1), ja(1), nloc, nrhs, bandwidth, iter, maxit
  integer  spikeparm(1), num_procs, rank
  double precision a(1), tol, B(nloc,1), X(nloc,1)
  double precision P(nloc,1), Y(nloc,1), W(nloc,1), R(nloc,1)
  double precision  xtmp(1), xtmp_loc(1), ddot, dnrm2, dsqrt, T_cgmatvec, T_cg_allred, Tstart

  integer i, j, job, ierr
  double precision stagtest, ddum 
 
  double precision, dimension(:), allocatable :: prn, prn0, alpha, beta, mu, rho, rho_old

!=============================================================
  allocate(prn(nrhs), prn0(nrhs), alpha(nrhs), beta(nrhs), mu(nrhs), rho(nrhs), rho_old(nrhs))


  do j = 1, nrhs
     do i = 1, nloc;  X(i,j) = 0.0d0; R(i,j) = B(i,j); enddo
  enddo
  do j = 1, nrhs
     do i = 1, nloc;  W(i,j) = R(i,j);   enddo
  enddo
  do j = 1, nrhs
     alpha(j) = 0.0d0;   beta(j) = 0.0d0;
  enddo

if (prec) then 
  ! Z =   M \ R;  
  job = 88; !! JOB = 2 for solve
  if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
  call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MATVEC(job, neqns, nzmax, ia, ja, a, W,nloc, nrhs, num_procs, rank )
end if 

!do j = 1, nrhs;
!     xtmp_loc(j) = dnrm2(nloc, W(1,j), 1)**2
!enddo
!call MPI_Allreduce(xtmp_loc, prn0, nrhs, MPI_DOUBLE_PRECISION,  &
!      MPI_SUM, MPI_COMM_WORLD, ierr)
do j = 1, nrhs;
     xtmp_loc(j) = maxval(abs(W(1:nloc,j))) 
enddo
call MPI_Allreduce(xtmp_loc, prn0, nrhs, MPI_DOUBLE_PRECISION,  &
      MPI_MAX, MPI_COMM_WORLD, ierr)





do iter = 1,maxit 
  ! rho = R' * Z ;
!$OMP PARALLEL DO PRIVATE(j)
  do j = 1, nrhs
     xtmp_loc(j) = ddot(nloc,R(1,j),1,W(1,j),1)
  enddo
!$OMP END PARALLEL DO
  Tstart = MPI_Wtime()
  call MPI_Allreduce(xtmp_loc, rho, nrhs, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)
  T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart

if (iter .eq. 1) then 
  ! P = Z ;
!$OMP PARALLEL DO PRIVATE(i,j)
  do j = 1, nrhs
   do i = 1, nloc;  P(i,j) = W(i,j);   enddo
  enddo
!$OMP END PARALLEL DO
else 
   !! beta = rho / rho_old ;
!$OMP PARALLEL DO PRIVATE(j)
   do j = 1, nrhs
     beta(j)  = rho(j)/rho_old(j);
   end do
!$OMP END PARALLEL DO
   ! p = z + beta*p ; 
!$OMP PARALLEL DO PRIVATE(i,j)
     do j = 1, nrhs
        do i = 1, nloc
            P(i,j) = W(i,j) + beta(j) * P(i,j);
        enddo
     enddo
!$OMP END PARALLEL DO
end if 

!$OMP PARALLEL DO PRIVATE(j)
do j=1,nrhs 
   rho_old(j) = rho(j) 
end do 
!$OMP END PARALLEL DO

  ! Ap = A*p;  
!$OMP PARALLEL DO PRIVATE(i,j)
     do j = 1, nrhs
        do i = 1, nloc
            Y(i,j) = P(i,j);
        enddo
     enddo
!$OMP END PARALLEL DO
  Tstart = MPI_Wtime()
  job = 55; !! JOB = 2 for solve
  if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
  call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MATVEC(job, neqns, nzmax, ia, ja, a, Y,nloc, nrhs, num_procs, rank )
  T_cgmatvec = T_cgmatvec + MPI_Wtime() - Tstart

  ! mu = P' * AP 
!$OMP PARALLEL DO PRIVATE(j)
  do j = 1, nrhs
     xtmp_loc(j) = ddot(nloc,P(1,j),1,Y(1,j),1)
  enddo
!$OMP END PARALLEL DO

  Tstart = MPI_Wtime()
  call MPI_Allreduce(xtmp_loc, mu, nrhs, MPI_DOUBLE_PRECISION,  &
                   MPI_SUM, MPI_COMM_WORLD, ierr)
  T_cg_allred = T_cg_allred + MPI_Wtime() - Tstart

  ! alpha = rho / mu 
!$OMP PARALLEL DO PRIVATE(j) 
  do j = 1, nrhs 
     alpha(j)  = rho(j)/mu(j); 
  end do 
!$OMP END PARALLEL DO




!$OMP PARALLEL DO PRIVATE(i,j)
     do j = 1, nrhs
        do i = 1, nloc
           X(i,j) = X(i,j) + alpha(j) * P(i,j)
           R(i,j) = R(i,j) - alpha(j) * Y(i,j)
           W(i,j) = R(i,j) 
         enddo
     enddo
!$OMP END PARALLEL DO


! Z = D\R; 
if (prec) then 
  job = 88; !! JOB = 2 for solve
  if (rank .eq. 0 ) nzmax = ia(neqns+1)-1
  call MPI_BCAST(nzmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MATVEC(job, neqns, nzmax, ia, ja, a, W,nloc, nrhs, num_procs, rank )
end if 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!Check the relative residual   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP PARALLEL DO PRIVATE(j)
do j = 1, nrhs;
     xtmp_loc(j) =  maxval(abs(W(1:nloc,j))) 
enddo
!$OMP END PARALLEL DO 
call MPI_Allreduce(xtmp_loc, prn, nrhs, MPI_DOUBLE_PRECISION,  &
      MPI_MAX, MPI_COMM_WORLD, ierr) 
!$OMP PARALLEL DO PRIVATE(j)
do j =1,nrhs
     prn(j) =  prn(j) / prn0(j) 
end do 
!$OMP END PARALLEL DO 



!if (rank.eq. 0 ) write(*,*) prn 
if (maxval(prn) .le. tol ) exit 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! if not converging in 10 iteration exit
!if ((iter .eq. 10) .and. (maxval(prn) .ge. 1.0 )) exit


end do 

if (rank .eq. 0 ) write(*,*) 'CG (max) relative residual:', iter, maxval(prn)

!     do j = 1, nrhs;
!        xtmp_loc(j) = dnrm2(nloc, R(1,j), 1)**2
!     enddo
!     call MPI_Allreduce(xtmp_loc, prn, nrhs, MPI_DOUBLE_PRECISION,  &
!                   MPI_SUM, MPI_COMM_WORLD, ierr)
!
!     if (rank .eq. 0) write(*,*) 'PCG residual: ', sqrt(prn) 


return
end subroutine pcg

