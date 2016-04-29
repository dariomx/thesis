!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the equivalent of general_fiedler() for the case of matrices
! with >= 1 strongly connected components in their corresponding graphs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mc_general_fiedler(n, nnz_in, Ai_in, Aj_in, Av_in, perm, &
               max_outer_it, max_inner_it, inner_tol, outer_tol, abs_res_tol, lam_tol, &
               rank, num_procs)

  use mpi 
  use omp_lib
  
  implicit none
  
  integer, intent(in) :: n, nnz_in
  integer, dimension(1:nnz_in), intent(in) :: Ai_in, Aj_in
  double precision, dimension(1:nnz_in), intent(in) :: Av_in

  integer, dimension(1:n), intent(out) :: perm
  integer, intent(in) :: rank, num_procs

  integer max_outer_it, max_inner_it
  double precision   inner_tol, outer_tol, abs_res_tol, lam_tol
  

  integer nnz
  integer, allocatable, dimension(:) :: Ai, Aj
  double precision, allocatable, dimension(:) :: Av

  integer, allocatable, dimension(:) :: iAcsr, jAcsr
  double precision, allocatable, dimension(:) :: Acsr

  integer, allocatable, dimension(:) :: tmp_iA, tmp_jA, tmp_iAT, tmp_jAT
  double precision, allocatable, dimension(:) :: tmp_A, tmp_AT

  integer, allocatable, dimension(:) :: iwork

  integer :: nelems
  integer :: nb
  integer, dimension(:), allocatable :: p
  integer, dimension(:), allocatable :: r
  
  integer, dimension(:), allocatable :: brows
  integer, dimension(:), allocatable :: bcols
  double precision, dimension(:), allocatable :: bvals
  integer, dimension(:), allocatable :: brangelimits
  integer, dimension(:), allocatable :: bnnzlimits
  integer, dimension(:), allocatable :: bsizes
  integer, dimension(:), allocatable :: Nj_list, pstack
  integer, dimension(:), allocatable :: Aij, Ajj
  double precision, dimension(:), allocatable :: Avj
  integer :: i1, i2, k1, k2, Nj, nnzj
  
  integer :: trivial_end, dense_end
  integer, dimension(:), allocatable :: permj, perminv

  integer, dimension(:), allocatable :: perminvj

  integer :: ierr
  integer :: i, j, k, kk
  double precision dnrm2

  
  if (rank .eq. 0) then
     !=======================
     !
     !   nnz_in, Ai_in, Aj_in, Av_in  --> nnz, Ai, Aj, Av
     !
       allocate(tmp_A(nnz_in), tmp_iA(nnz_in), tmp_jA(nnz_in)) 
       allocate(tmp_AT(nnz_in), tmp_jAT(nnz_in), tmp_iAT(nnz_in)) 

       call coocsr(n,nnz_in, Av_in, Ai_in, Aj_in,  tmp_A, tmp_jA, tmp_iA)
       call coocsr(n,nnz_in, Av_in, Aj_in, Ai_in,  tmp_AT, tmp_jAT, tmp_iAT)

       allocate(iwork(max(n +1 , 2*nnz_in)))
       call csort(n, tmp_A, tmp_jA, tmp_iA, iwork, .true. )
       call csort(n, tmp_AT, tmp_jAT, tmp_iAT, iwork, .true. )
       deallocate(iwork)
       
       tmp_A = abs(tmp_A);        tmp_AT = abs(tmp_AT) 

       allocate(Acsr(2*nnz_in), iAcsr(n+1), jAcsr(2*nnz_in))

       call aplb1(n,n,1, tmp_A,tmp_jA,tmp_iA, tmp_AT, tmp_jAT, tmp_iAT,  &
                          Acsr, jAcsr, iAcsr, 2*nnz_in,ierr)

       deallocate(tmp_A, tmp_iA, tmp_jA)
       deallocate(tmp_AT,tmp_iAT,tmp_jAT)

       if (ierr .ne. 0 ) then
          write(*,*) 'aplb1 error:', ierr;     stop
       end if
       Acsr = Acsr/2.0

       nnz = iAcsr(n+1) - 1
       allocate(Av(nnz), Ai(nnz), Aj(nnz))

       kk = 0
       do i = 1, n
          do k = iAcsr(i), iAcsr(i+1)-1
             j = jAcsr(k)
             kk = kk+1
             Av(kk) = Acsr(k)
             Ai(kk) = i;        Aj(kk) = j
          enddo
       enddo

       print *, ' nnz, kk should be equal:  ', nnz, kk

       deallocate(Acsr, iAcsr, jAcsr)
     !
     !  Work with the newly created COO matrix Av, Ai, Aj (and nnz)
     !  Deallocate(Av, Ai, Aj) at the end of this routine
     !  
     !=======================

     nelems = nnz
     allocate(p(n))
     allocate(r(n+1))
     allocate(brows(nelems))
     allocate(bcols(nelems))
     allocate(bvals(nelems))
     allocate(brangelimits(n+1))
     allocate(bnnzlimits(n+1))
     
     call sccblocks(Ai, Aj, Av, n, nelems, nb, p, r, brows, bcols, bvals, brangelimits, bnnzlimits)
     call method_limits(nb, brangelimits, bnnzlimits, trivial_end, dense_end)

     print *, ' Number of components, N = ', nb, n
     
     allocate(Nj_list(nb))
     do j = 1, nb
        Nj_list(j) = brangelimits(j+1) - brangelimits(j)
     enddo
     allocate(pstack(n))
     allocate(perminv(n))
  end if
  
  call MPI_BCAST(nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(trivial_end, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(dense_end, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
  if (rank .ne. 0) allocate(Nj_list(nb))

  call MPI_BCAST(Nj_list, nb, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  
  do j = 1, nb
     Nj = Nj_list(j)

     if ((rank .eq. 0)) then
        write(*,"(a,i10,a,i10,a)") 'Component ', j, ' has ', Nj, ' vertices.'

        allocate(permj(Nj))
        nnzj = bnnzlimits(j+1) - bnnzlimits(j)
        allocate(Avj(nnzj), Aij(nnzj), Ajj(nnzj))
        k1 = bnnzlimits(j)
        k2 = bnnzlimits(j+1)-1
        Avj(1:nnzj) = bvals(k1:k2)
        Aij(1:nnzj) = brows(k1:k2)
        Ajj(1:nnzj) = bcols(k1:k2)
     end if

     if ((Nj .le. trivial_end) .and. (rank .eq. 0)) then
        call general_fiedler_trivial(Nj, nnzj, Aij, Ajj, Avj, permj)
     else if ((Nj .le. dense_end) .and. (rank .eq. 0)) then
        call general_fiedler_dense(Nj, nnzj, Aij, Ajj, Avj, permj)
     else if (Nj .gt. dense_end) then
        call general_fiedler(Nj, nnzj, Aij, Ajj, Avj, permj, &
           max_outer_it, max_inner_it, inner_tol, outer_tol, abs_res_tol, lam_tol, &
           rank, num_procs)
     end if
     
     if (rank .eq. 0) then

        allocate(perminvj(Nj))

        do i = 1, Nj
           perminvj(permj(i)) = i
        enddo

        i1 = brangelimits(j)
        i2 = brangelimits(j+1)-1

!!!        pstack(i1:i2) = permj(1:Nj)
        pstack(i1:i2) = perminvj(1:Nj)


        deallocate(Avj, Aij, Ajj, permj)
        deallocate(perminvj)
     end if
  enddo
  
  if (rank .eq. 0) then
     call calcperm(p, pstack, n, nb, brangelimits, perminv)

     !  Alicia's fix
     do i = 1, n
        perm(perminv(i)) = i
     enddo
  end if


!    conditional added to debug:
  if (rank .eq. 0) then
     deallocate(Av, Ai, Aj)
  endif


return

end subroutine mc_general_fiedler






