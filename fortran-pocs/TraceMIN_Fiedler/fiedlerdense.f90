!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the matrix size(order) limits for the application of the three permutation computing
! methods in mc_general_fiedler():
! If size <= trivial_end then use general_fiedler_trivial()
! If trivial_end < size <= dense_end then use general_fiedler_dense()
! If size > dense_end use then use general_fiedler() (i.e. the (parallel) TracerMin-Fiedler routine). 
!
! XXXX A more educated guess could be made for the limits if the distributions in the nb blocks 
! are taken into account, i.e. the block matrix sizes that can be accessed from brangelimits and
! the number of nonzero elements of each block that can be accessed from bnnzlimits.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine method_limits(nb, brangelimits, bnnzlimits, trivial_end, dense_end)
  implicit none
  integer, intent(in) :: nb
  integer, intent(in) :: brangelimits(nb+1), bnnzlimits(nb+1)
  integer, intent(out) :: trivial_end, dense_end
  
  trivial_end = 6   !!  2 !! 2
  dense_end =  80  !!  200  !!  500
end subroutine method_limits



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the trivial identity permutation for a matrix geiven in coordinate format.
! It does not use any Fiedler machinery or matrix information except for its order n.
! Redundant argument are there just for the uniformity of calling conventions
! for the three permutation computing methods currently available.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine general_fiedler_trivial(n, nnz, Ai, Aj, Av, perm)
  implicit none
  integer, intent(in) :: n, nnz
  integer, intent(in) :: Ai(nnz), Aj(nnz)
  double precision, intent(in) :: Av(nnz)
  integer, intent(out) :: perm(n)
  
  integer :: i
  
  do i=1,n
     perm(i) = i
  end do

end subroutine general_fiedler_trivial



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the permutation perm induced by the Fiedler vector as computed by the Laplacian of a
! matrix given in coordinate format.
! Internally the matrix is first converted to dense format, then its Laplacian is computed
! and finally the Fiedler vector and the permutation are generated.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine general_fiedler_dense(n, nnz, Ai, Aj, Av, perm) 
  implicit none
  integer, intent(in) :: n, nnz
  integer, intent(in) :: Ai(nnz), Aj(nnz)
  double precision, intent(in) :: Av(nnz)
  integer, intent(out) :: perm(n)
  
  double precision, dimension(:,:), allocatable :: A
  double precision, dimension(:,:), allocatable :: L
  integer :: k, kmax = 6
  double precision, dimension(:), allocatable :: eigenvals
  double precision, dimension(:,:), allocatable :: eigenvecs
  integer :: error
  double precision, dimension(:), allocatable :: fiedlervec
  integer, dimension(:), allocatable :: perminv
  integer :: i
  
  ! Compute Laplacian
  allocate(A(n, n))
  allocate(L(n, n))
  
  call cootodense(Ai, Aj, Av, n, nnz, A)
  call laplacian_dense(A, n, L)
  
  ! Compute k eigenpairs
  k = min(kmax, n)
  allocate(eigenvals(k))
  allocate(eigenvecs(n, k))
  call eigen(L, n, k, eigenvals, eigenvecs, error)
  
  ! Compute Fiedler-induced permutation
  allocate(perminv(n))
  allocate(fiedlervec(n))
  fiedlervec(1:n) = eigenvecs(1:n, 2)
  call idxsort(fiedlervec, n, perminv)
  do i=1,n
     perm(perminv(i)) = i
  end do
  
  ! Clean-up
  deallocate(fiedlervec, perminv, eigenvecs, eigenvals, L, A)
  
end subroutine general_fiedler_dense



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert a matrix input in coordinate format into a dense matrix format returned in A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cootodense(rows, cols, vals, n, nnz, A)
  implicit none
  integer, intent(in) :: rows(nnz), cols(nnz)
  double precision, intent(in) :: vals(nnz)
  integer, intent(in) :: n, nnz
  double precision, intent(out) :: A(n, n)
  
  integer :: i
  A = 0.0
  do i=1,nnz
     A(rows(i), cols(i)) = vals(i)
  end do
end subroutine cootodense



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the Laplacian of square dense matrix A of order n; the result will be in L
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine laplacian_dense(A, n, L)
  implicit none
  double precision, intent(in) :: A(n, n)
  integer, intent(in) :: n
  double precision, intent(out) :: L(n, n)
  
  integer :: i, j
  
  L = -(A + transpose(A))
  ! where (L .ne. 0.0 ) L = -1.0
  do i=1,n
     ! L(i, i) = count(L(i,1:n) .ne. 0.)
     L(i, i) = abs(sum(L(i,1:n))) - abs(L(i,i))
  end do
end subroutine laplacian_dense



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the first k eigenvectors of square matrix L of order n 
! The eignevalues are put in eigenvals (array of size k) and the eigenvectors
! appear as columns of eigenvecs (array of size n x k).
! error copies the value of info variable in lapack routine dsyevx() 
! (that is actually used). 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eigen(L, n, k, eigenvals, eigenvecs, error)
  implicit none
  double precision, intent(in) :: L(n, n)
  integer, intent(in) :: n
  integer, intent(in) :: k
  double precision, intent(out) :: eigenvals(k)
  double precision, intent(out) :: eigenvecs(n, k)
  integer, intent(out) :: error
  
  character :: jobz, range, uplo 
  integer :: lda, il, iu, m, ldz, lwork, info
  double precision :: vl , vu, abstol
  double precision, dimension(:), allocatable :: w, work
  double precision, dimension(:, :), allocatable :: Z
  integer, dimension(:), allocatable :: iwork, ifail
  
  jobz = 'V'
  range = 'I'
  uplo = 'U'
  lda = n
  vl = 0.0
  vu = 0.0
  il = 1
  iu = k
  abstol = 1.0e-7
  ldz = lda
  lwork = 8 * n
  allocate(w(n))
  allocate(Z(ldz, iu-il+1))
  allocate(work(lwork))
  allocate(iwork(5*n))
  allocate(ifail(n))
  call dsyevx(jobz, range, uplo, n, L, lda, vl, vu, il, iu, abstol, &
       m, w, Z, ldz, work, lwork, iwork, ifail, info)
  if (info .eq. 0) then
     eigenvals(1:k) = w(1:k)
     eigenvecs(1:n,1:k) = Z(1:n,1:k)
  end if
  error = info
  
  ! Clean-up
  deallocate(w, work, Z, iwork, ifail)

end subroutine eigen



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the permutation that links to the original graph node labels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calcperm(p, pstack, n, nb, brangelimits, perminv)
  implicit none 
  integer, intent(in) :: p(n), pstack(n), brangelimits(nb+1)
  integer, intent(in) :: n, nb
  integer, intent(out) :: perminv(n)
  
  integer :: i, j, k, bstart, bend,   ii
  integer, dimension(:), allocatable :: pnew

  allocate(pnew(n))
  
  do i=1,nb
     bstart = brangelimits(i)
     bend = brangelimits(i+1)
     do j=bstart, bend-1
        k = bstart + pstack(j) - 1
        pnew(j) = p(k)
     end do
  end do

  do i = 1, n
!     perminv(pnew(i))=i
     perminv(i)=pnew(i)
  enddo

! debug print
!  do i = 1, n
!     write(*,"(4i10)") i, p(i), pstack(i), perminv(i)
!  enddo

  
  ! Clean-up
  deallocate(pnew)

end subroutine calcperm



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write a matrix in coordinate format to a file in matrix market format (.mtx)
! (coordinate, real, general)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mtxwrite_matrix(name, n, nnz, Ai, Aj, Av)
  implicit none
  
  character name*80 
  integer, intent(in) :: n, nnz
  integer, intent(in) :: Ai(nnz), Aj(nnz)
  double precision, intent(in) :: Av(nnz)
  
  integer ounit
  character rep*10
  character field*7
  character symm*19
  integer, dimension(:), allocatable :: IDummy
  double complex, dimension(:), allocatable :: CDummy
  
  
  ounit = 24
  rep = 'coordinate'
  field = 'real'
  symm = 'general'
  
  allocate(IDummy(nnz))
  allocate(CDummy(nnz))
  
  open(unit=ounit, file=name)
  call mmwrite(ounit, rep, field, symm, n, n, nnz, Ai, Aj, IDummy, Av, CDummy)
  close(ounit)

end subroutine mtxwrite_matrix



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write a vector of integers in coordinate format to a file in matrix market format (.mtx)
! (coordinate, integer, general), write it as a single column matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mtxwrite_vector(name, n, v)
  implicit none
  
  character name*80 
  integer, intent(in) :: n
  double precision, intent(in) :: v(n)
  
  integer ounit
  character rep*10
  character field*7
  character symm*19
  double precision nnz
  integer, dimension(:), allocatable :: vi, vj
  double precision, dimension(:), allocatable :: RDummy
  double complex, dimension(:), allocatable :: CDummy
  integer :: i
  

  
  ounit = 25
  rep = 'coordinate'
  field = 'integer'
  symm = 'general'
  
  allocate(RDummy(n))
  allocate(CDummy(n))
  allocate(vi(n))
  allocate(vj(n))
  vj = 1
  do i=1,n
     vi(i) = i
  end do
  
  open(unit=ounit, file=name)
  call mmwrite(ounit, rep, field, symm, n, 1, n, vi, vj, v, RDummy, CDummy)
  close(ounit)

end subroutine mtxwrite_vector


