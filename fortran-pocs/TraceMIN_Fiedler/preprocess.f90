!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Computes the Fiedler vector using the TraceMin algorithm                !
!                                                                            !
!    Copyright (C) 2010 Murat Manguoglu, Faisal Saied, Eric Cox, Ahmed Sameh,!
!    Purdue University                                                       !
!    Contact: murat.manguoglu@gmail.com                                      !
!                                                                            !
!    This file is a part of TraceMin-Fiedler                                 !
!    TraceMin-Fiedler is free software: you can redistribute it and/or modify!
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
!                                                                            !
!    By downloading the software TraceMin-Fiedler I agree to give reference  !
!    to the authors of the software in any publications resulting from its   !
!    use.                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine general_fiedler(n , nnz , Ai , Aj, Av,perm, &
      max_outer_it, max_inner_it, inner_tol, outer_tol, abs_res_tol, lam_tol, &
      rank,num_procs) 

        use mpi 
        use omp_lib
        implicit none


       integer ::ierr, error,  k,  i,j , nnz,nz , n ,m, cnt, rank, num_procs, max_outer_it , max_inner_it  
       double precision :: t0, t1 , MC73_TOL, MC73_TOL1, MC73_RTOL, total , di
       double precision :: inner_tol, outer_tol, abs_res_tol, lam_tol

       double precision , dimension(:),allocatable ::ATA, Apv ,A , Ap, AT, Bv ,B 
       double precision,   dimension(:), allocatable :: fiedler_vec
       double precision  :: Av(nnz) 
       integer :: Ai(nnz), Aj(nnz) , perm(n) 
       integer, dimension(:),allocatable ::new_ptr , tmp, temp, iATA, jATA ,jAT, iAT, Bi, Bj , iB , jB, Api, Apj, iA,jA,iAp, jAp, perm_inv
       character(len=100) :: name

       logical, dimension(:),allocatable ::  diagonal

        integer rp(200)


      character rep*10
      character field*7
      character symm*19
      character tmp1*1024
      character tmp2*2
      integer,dimension(:),allocatable   :: IDummy
      double complex, dimension(:),allocatable   :: CDummy


      integer :: ICNTL(10) , NUM , LIW, LDW
      integer , dimension(:), allocatable:: IWORK, LENC
      double precision, dimension(:), allocatable:: DWORK
      double precision  CNTL(5)

      logical remove_dangling_nodes

if (rank .eq. 0) then 

       allocate(A(nnz)) 
       allocate(AT(nnz))
       allocate(jA(nnz))
       allocate(jAT(nnz))
       allocate(iA(n+1))
       allocate(iAT(n+1))

        call coocsr(n,nnz, Av, Ai, Aj, A,jA,iA)
        call coocsr(n,nnz, Av, Aj, Ai, AT,jAT,iAT)



       allocate(iwork(max(n +1 , 2*nnz)))
       call csort(n, A, jA, iA, iwork, .true. )
       deallocate(iwork)


       allocate(iwork(max(n +1 , 2*nnz)))
       call csort(n, AT, jAT, iAT, iwork, .true. )
       deallocate(iwork)

       A = abs(A)
       AT = abs(AT) 

        allocate(ATA(2*nnz))
        allocate(iATA(n+1))
        allocate(jATA(2*nnz))

!!!!!!!!!!!!!!!!!!!!!!!!!!
      call aplb1(n,n,1,A,jA,iA,AT,jAT,iAT,ATA,jATA,iATA,2*nnz,ierr)
      if (ierr .ne. 0 ) then
          write(*,*) 'aplb1 error:', ierr
         stop
      end if
      ATA = ATA/2.0
!!!!!!!!!!!!!!!!!!!

      allocate(diagonal(n))
      allocate(new_ptr(n))
    
      diagonal =  .false.

      remove_dangling_nodes = .false.
!!!      remove_dangling_nodes = .true.

      if (remove_dangling_nodes) then
         print *, ' general_fiedler in preprocess.f90:  Dangling nodes removed'
         cnt = 0
         do i =1, n
            total = 0.0d0
            di = 0.0d0
            do  k=iATA(i) , iATA(i+1)-1
               total = total + abs(ATA(k))
               if (i .eq. jATA(k))  di = abs(ATA(k))
            end do
            if (total .eq. di) then
               diagonal(i) = .true.
               ! write(*,*) i,i, ' is diagonal' 
               cnt = cnt + 1
            end if
            
         end do
      else
         print *, ' general_fiedler in preprocess.f90:  Dangling nodes not removed'
         print *, ' general_fiedler in preprocess.f90:  N must equal cnt '
      endif
         
      new_ptr = 0
      cnt = 0
      do i =1,n
      if (.not. diagonal(i)) then
        cnt = cnt +1
        new_ptr(i) = cnt
      end if
      end do

       nz = 0
       do i=1, n
          do  k=iATA(i) , iATA(i+1)-1
              if (i .ge. jATA(k)) then
                if (.not. diagonal(jATA(k)) .and. .not. diagonal(i)) nz = nz +1
              end if
          end do
       end do
deallocate(A,iA,jA)
deallocate(AT,iAT, jAT)

allocate(jB(nz)) 
allocate(B(nz)) 
allocate(iB(cnt+1)) 

allocate(Bi(nz))
allocate(Bj(nz))
allocate(Bv(nz))

nz = 0 
do i=1, n 
   do  k=iATA(i) , iATA(i+1)-1 
      if (i .ge. jATA(k)) then 
         if (.not. diagonal(jATA(k)) .and. .not. diagonal(i)) then
            nz = nz + 1 
            Bi(nz) = new_ptr(i) 
            Bj(nz) = new_ptr(jATA(k))
            Bv(nz) = ATA(k)
         end if
      end if
   end do
end do

call coocsr(cnt,nz, Bv, Bi, Bj, B,jB,iB)
deallocate(Bv, Bi, Bj) 
deallocate(ATA,jATA,iATA)

end if

call MPI_BCAST(cnt , 1, MPI_INTEGER,0,MPI_COMM_WORLD,error)


if (rank .eq. 0) allocate(fiedler_vec(cnt))


!  These are now passed in and not defined here
!max_outer_it = 150
!max_inner_it = 6000
!outer_tol = 1.0d-9  !!  1.0d-10   !!  1.0d-7
!inner_tol = 1.0d-11   !! 4
!
! Also, abs_res_tol, lam_tol are passed in.
!

if (rank .eq. 0) then
   print *, ' Calling Fiedler:  N, cnt = ', N, cnt
endif


call FIEDLER(cnt, nz, iB, jB, B,fiedler_vec,max_outer_it, max_inner_it , outer_tol, inner_tol,  &
      abs_res_tol, lam_tol, num_procs, rank )


if (rank .eq. 0 ) then
deallocate(iB,jB,B) 
allocate(tmp(cnt)) 
allocate(perm_inv(n)) 

!------------------------------------------------
!
! replace kb07ad by idxsort
!
!!!!!! call KB07AD (fiedler_vec, cnt, tmp) 

call idxsort(fiedler_vec, cnt, tmp) 
!
!------------------------------------------------
 do i = 1, n 
   if (.not. diagonal(i)) then   
    perm_inv(i) = tmp(new_ptr(i))
   else 
     cnt = cnt + 1 
     perm_inv(i) = cnt 
   end if 
 end do

  do i = 1, n
    perm(perm_inv(i)) = i
  end do
deallocate(perm_inv) 
deallocate(tmp)
deallocate(diagonal) 
deallocate(new_ptr) 
deallocate(fiedler_vec) 

end if

end subroutine general_fiedler 

