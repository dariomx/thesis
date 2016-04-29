!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Computes the Laplacian, distributes it, performs matrix-vector          !
!    multiplications, diagonal precondtioning                                !
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
!  									     !
!    By downloading the software TraceMin-Fiedler I agree to give reference  !
!    to the authors of the software in any publications resulting from its   !
!    use.							             ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine MATVEC(step, neqns, nzmax, ia, ja, a, f,nloc, nrhs, nb_procs, rank )
        use mpi
        use omp_lib
	implicit none 

	! Function parameter declarations
	integer step
        integer nloc 
        integer ia(neqns+1) , ja(nzmax) 
        double precision a(nzmax) 
	integer  kk, neqns, nzmax, nrhs 
	double precision f(nloc,nrhs)
	double precision tol, nrm 
	double precision, save :: perturb,  bicg_time = 0
	integer  k, j,  nb_procs, rank , spikeparm(100) 

	! Static variables
        integer :: nrequest(nb_procs), request 
        integer, dimension(:), allocatable :: iB, jB 
        double precision, dimension(:), allocatable :: B 
	integer, save :: i1, i2, j1, j2, ratio, kmax, nr, nc, nz, max_it, msglvl, maxfct, mtype, mnum,n,  phase
        integer, save :: part1, part1_LAST 
	integer*8, save :: pt(64)
	integer, dimension(:), allocatable, save ::iA2, jA2, Ai, Aj, rcounts, displs, iBBB, jBBB, iAloc, jAloc, prec_iA, prec_jA, permp,permu, IPIV2, iBB, jBB, iBBTBB, jBBTBB, iAA, jAA
	integer, dimension(:), allocatable, save :: perm,perm_inv,  bja, bia, iwork , IIWORK 
	double precision, dimension(:), allocatable, save :: A2, Av, BB, BBB, BBTBB, Aloc, x1, x2, f1, f2, u1, u2, D,Di,  prec_A, Dhalf, AA
	double precision, dimension(:,:), allocatable, save :: BC, RED, VW, SENDR, RECVR, BND, SENDK, RECVK, p, p_hat, v, s_hat, r, s, t, xi, r_old, SEND1, SEND2, RECV1, RECV2 
	double precision, dimension(:), allocatable, save :: scaled , DDWORK, scaling, ba
        double precision, dimension(:,:), allocatable, save :: ftmp ,temp1, temp2, temp, tmp  
        integer, dimension(:), allocatable, save :: perm_mc64, perm_mc64_inv  
	integer, dimension(:), save :: iparm(64)
        integer, save:: count=0
        logical, save :: mc64

!!! MKL tips related
         integer, save :: uborder , lborder
!!!!!!!!!!!!!!!!!!!!!
!!! ILUT variables
      double precision,save, dimension(:), allocatable :: alu, wu , wl
      integer ,save,  dimension(:), allocatable:: iw,  jlu , ju , jr, jwu , jwl
      integer iwk, LFILL
      double precision DROPTOL
!!!!!!!!!!!!!!!!!


        double precision ::norm(nrhs),  bnrm2(nrhs) , omega(nrhs) , omega1(nrhs) , alpha(nrhs) , beta(nrhs), rho(nrhs), relres(nrhs), rho_1(nrhs) 


	! Local variables
	double precision t0 , t1
	integer  i, k2, code, error, alloc_error
        double precision, dimension(:,:), allocatable ::rhs
 
       
	!!!!!!!!!!! unused local variables 
	integer  prec_nnz, nnz,  l, prec_k_times, kl, ku, cnt, req(4), idum, status2(MPI_STATUS_SIZE, 4), LDW, LIW , NUM
         
	integer, dimension(MPI_STATUS_SIZE) :: status
        integer ICNTL(10), INFO(10) 
	double precision ddot, ddum(nrhs), RINFO(20)

	integer :: partition, cut_size = 100 



logical, dimension(:,:) , allocatable,save :: map
integer, dimension(:), allocatable,save ::  n_send
integer, dimension(:) , allocatable,save :: n_recv
double precision, dimension(:) ,allocatable,save :: SENDBUF, RECVBUF
integer, dimension(:), allocatable,save :: reqsend, reqrecv
logical, dimension(:), allocatable,save :: recvd
integer , dimension(:,:) , allocatable,save :: sendmap, recvmap



if (step .eq. 88) then  ! matvec + diagonal preconditioner  

!$OMP PARALLEL DO PRIVATE(i)
               do i = 1, nr
                      f(i,1:nrhs) = f(i,1:nrhs)/Di(i) 
               end do
!$OMP END PARALLEL DO

end if 

if (step .eq. 77) then  ! matvec + diagonal preconditioner

!$OMP PARALLEL DO PRIVATE(i)
               do i = 1, nr
                      f(i,1:nrhs) = f(i,1:nrhs)/(Di(i)+perturb)
               end do
!$OMP END PARALLEL DO

end if



if (step .eq. 55) then 
call dcopy(nr*nrhs,f,1,p_hat(i1:i2,1:nrhs),1) 
    do i=1, nb_procs
      if (n_send(i) .ne. 0) then
           do j=1,n_send(i)
              do k=1,nrhs
                SENDBUF(k+(j-1)*nrhs+(i-1)*nrhs*maxval(n_send)) = p_hat(sendmap(j,i),k)
              end do
           enddo
           call MPI_ISEND(SENDBUF(1+(i-1)*nrhs*maxval(n_send)),nrhs*n_send(i),MPI_DOUBLE_PRECISION,i-1,555+rank+1,MPI_COMM_WORLD,reqsend(i),error)
      end if
      if(n_recv(i) .ne. 0 ) then
          call MPI_IRECV(RECVBUF(1+(i-1)*nrhs*maxval(n_recv)),nrhs*n_recv(i), MPI_DOUBLE_PRECISION, i-1,555+i, MPI_COMM_WORLD,reqrecv(i) ,error)
      end if
    end do

    do i =1, nb_procs
      if ( n_recv(i).ne. 0 ) then
         call MPI_WAIT(reqrecv(i),status,error)
         do j=1,n_recv(i)
            do k=1,nrhs
             p_hat(recvmap(j,i),k) = RECVBUF(k+(j-1)*nrhs+(i-1)*nrhs*maxval(n_recv))
            end do
         enddo
      end if
    end do

!$OMP PARALLEL DO PRIVATE(i)
   do i = 1, nr
        f(i,1:nrhs) = 0.0d0
        do  k = iAloc(i), iAloc(i+1)-1
           f(i,1:nrhs) = f(i,1:nrhs) +  Aloc(k) * p_hat(jAloc(k),1:nrhs)
        end do
   end do
!$OMP END PARALLEL DO


end if 




if (step .eq. 44) then
call dcopy(nr*nrhs,f,1,p_hat(i1:i2,1:nrhs),1)

    do i=1, nb_procs
      if (n_send(i) .ne. 0) then
           do j=1,n_send(i)
              do k=1,nrhs
                SENDBUF(k+(j-1)*nrhs+(i-1)*nrhs*maxval(n_send)) = p_hat(sendmap(j,i),k)
              end do
           enddo
           call MPI_ISEND(SENDBUF(1+(i-1)*nrhs*maxval(n_send)),nrhs*n_send(i),MPI_DOUBLE_PRECISION,i-1,555+rank+1,MPI_COMM_WORLD,reqsend(i),error)
      end if
      if(n_recv(i) .ne. 0 ) then
          call MPI_IRECV(RECVBUF(1+(i-1)*nrhs*maxval(n_recv)),nrhs*n_recv(i), MPI_DOUBLE_PRECISION, i-1,555+i, MPI_COMM_WORLD,reqrecv(i) ,error)
      end if
    end do

    do i =1, nb_procs
      if ( n_recv(i).ne. 0 ) then
         call MPI_WAIT(reqrecv(i),status,error)
         do j=1,n_recv(i)
            do k=1,nrhs
             p_hat(recvmap(j,i),k) = RECVBUF(k+(j-1)*nrhs+(i-1)*nrhs*maxval(n_recv))
            end do
         enddo
      end if
    end do

 !!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP PARALLEL DO PRIVATE(i,k,k2)
   do i = 1, nr
        !f(i,1:nrhs) = 0.0d0
         f(i,1:nrhs) =  perturb  *p_hat(i,1:nrhs)
        do  k = iAloc(i), iAloc(i+1)-1
!!!!!!           f(i,1:nrhs) = f(i,1:nrhs) +  Aloc(k) * p_hat(jAloc(k),1:nrhs)
           do k2 = 1, nrhs
              f(i,k2) = f(i,k2) +  Aloc(k) * p_hat(jAloc(k),k2)
           enddo

        end do
   end do
!$OMP END PARALLEL DO

end if



if (step .eq. 66 )  then   !! sparse matvec
call dcopy(nr*nrhs,f,1,p_hat(i1:i2,1:nrhs),1)
if (ratio .ge. kmax) then
                if (rank .ne. nb_procs-1) then
                         call dcopy(nrhs*kmax,p_hat(1+((rank+1)*part1)-kmax:((rank+1)*part1),1:nrhs), 1, SEND1,1)
                         call MPI_ISEND(SEND1,nrhs*kmax,MPI_DOUBLE_PRECISION,rank+1,555+rank+1,MPI_COMM_WORLD,req(1),error)
                         call MPI_IRECV(RECV1,nrhs*kmax,MPI_DOUBLE_PRECISION,rank+1,999+rank+1, MPI_COMM_WORLD,req(3),error)
                end if
                if (rank .ne. 0 ) then
                            call dcopy(nrhs*kmax, p_hat(1+(rank)*part1:rank*part1+kmax,1:nrhs), 1, SEND2,1)
                            call MPI_ISEND(SEND2,nrhs*kmax,MPI_DOUBLE_PRECISION,rank-1,999+rank,MPI_COMM_WORLD,req(2), error)
                            call MPI_IRECV(RECV2,nrhs*kmax,MPI_DOUBLE_PRECISION,rank-1,555+rank, MPI_COMM_WORLD,req(4) ,error)
                end if
                if (rank .ne. 0) then
                    call MPI_WAIT(req(4), status,error)
                    call dcopy(kmax*nrhs,RECV2, 1, p_hat(1+((rank)*part1)-kmax:(rank)*part1,1:nrhs), 1)
                end if
                if(rank .ne. nb_procs-1) then
                  call MPI_WAIT(req(3), status,error)
                  call dcopy(kmax*nrhs, RECV1,1 , p_hat(1+(rank+1)*part1:(rank+1)*part1+kmax,1:nrhs) ,1 )
                end if


!$OMP PARALLEL DO PRIVATE(i)
               do i = 1, nr
                      f(i,1:nrhs) =  perturb  *p_hat(i,1:nrhs) 
                      !f(i,1:nrhs) = 0.0d0
                      do k = iAloc(i), iAloc(i+1)-1
                         f(i,1:nrhs) = f(i,1:nrhs) +  Aloc(k) *p_hat(jAloc(k),1:nrhs)
                      end do
                      !f(i,1:nrhs) = f(i,1:nrhs) + perturb  *p_hat(i,1:nrhs) 
               end do
!$OMP END PARALLEL DO
else
               do i=1,nrhs
               call MPI_ALLGATHERV(p_hat(i1:i2,i), nr ,MPI_DOUBLE_PRECISION, temp(:,i), rcounts, displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,error)
               end do
!$OMP PARALLEL DO PRIVATE(i)
               do i = 1, nr
                       !f(i,1:nrhs) = 0.0d0
                       f(i,1:nrhs) =  perturb  *temp(i,1:nrhs)
                       do k = iAloc(i), iAloc(i+1)-1
                         f(i,1:nrhs) = f(i,1:nrhs) +  Aloc(k) * temp(jAloc(k),1:nrhs)
                      end do
                      ! f(i,1:nrhs) = f(i,1:nrhs) + perturb  *temp(i,1:nrhs)
               end do
!$OMP END PARALLEL DO
end if
end if





if (step .eq. 99 )  then   !! sparse matvec 
call dcopy(nr*nrhs,f,1,p_hat(i1:i2,1:nrhs),1)
if (ratio .ge. kmax) then
                if (rank .ne. nb_procs-1) then
                         call dcopy(nrhs*kmax,p_hat(1+((rank+1)*part1)-kmax:((rank+1)*part1),1:nrhs), 1, SEND1,1)
                         call MPI_ISEND(SEND1,nrhs*kmax,MPI_DOUBLE_PRECISION,rank+1,555+rank+1,MPI_COMM_WORLD,req(1),error)
                         call MPI_IRECV(RECV1,nrhs*kmax,MPI_DOUBLE_PRECISION,rank+1,999+rank+1, MPI_COMM_WORLD,req(3),error)
                end if
                if (rank .ne. 0 ) then
                            call dcopy(nrhs*kmax, p_hat(1+(rank)*part1:rank*part1+kmax,1:nrhs), 1, SEND2,1)
                            call MPI_ISEND(SEND2,nrhs*kmax,MPI_DOUBLE_PRECISION,rank-1,999+rank,MPI_COMM_WORLD,req(2), error)
                            call MPI_IRECV(RECV2,nrhs*kmax,MPI_DOUBLE_PRECISION,rank-1,555+rank, MPI_COMM_WORLD,req(4) ,error)
                end if
                if (rank .ne. 0) then
                    call MPI_WAIT(req(4), status,error)
                    call dcopy(kmax*nrhs,RECV2, 1, p_hat(1+((rank)*part1)-kmax:(rank)*part1,1:nrhs), 1)
                end if
                if(rank .ne. nb_procs-1) then
                  call MPI_WAIT(req(3), status,error)
                  call dcopy(kmax*nrhs, RECV1,1 , p_hat(1+(rank+1)*part1:(rank+1)*part1+kmax,1:nrhs) ,1 )
                end if


!$OMP PARALLEL DO PRIVATE(i)
               do i = 1, nr
                      f(i,1:nrhs) = 0.0d0
                      do k = iAloc(i), iAloc(i+1)-1
                         f(i,1:nrhs) = f(i,1:nrhs) +  Aloc(k) *p_hat(jAloc(k),1:nrhs)
                      end do
               end do
!$OMP END PARALLEL DO
else
               do i=1,nrhs
               call MPI_ALLGATHERV(p_hat(i1:i2,i), nr ,MPI_DOUBLE_PRECISION, temp(:,i), rcounts, displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,error)
               end do
!$OMP PARALLEL DO PRIVATE(i)
               do i = 1, nr
                      f(i,1:nrhs) = 0.0d0
                       do k = iAloc(i), iAloc(i+1)-1
                         f(i,1:nrhs) = f(i,1:nrhs) +  Aloc(k) * temp(jAloc(k),1:nrhs)
                      end do
               end do
!$OMP END PARALLEL DO
end if


end if 


if ((step .eq. 1).or. (step .eq. 2))  then
if (rank .eq. 0 ) then 


!	allocate(perm(neqns)) 
!	allocate(perm_inv(neqns)) 
!do i=1, neqns 
!   perm(i) = i 
!   perm_inv(i) = i 
!end do 
!	call genrcm (neqns, nzmax, ia, ja, perm )
!	call perm_inverse3 (neqns, perm, perm_inv )
!	allocate(B(nzmax)) 
!	allocate(iB(neqns+1))  
!	allocate(jB(nzmax))  
!	call dperm(neqns, a, ja, ia, B, jB, iB,perm_inv,perm_inv, 3)
!        deallocate(perm_inv) 
!        deallocate(perm) 

      nnz = 0
      do i=1, neqns
          do  k=ia(i) , ia(i+1)-1
                 if (i .ne. ja(k)) then
                    nnz = nnz + 1
                 end if
          end do
      end do


      do i=1, neqns
          do  k=ia(i) , ia(i+1)-1
              if (i .ne. ja(k)) then
                  nnz = nnz + 1
              end if
          end do
      end do

      nnz = nnz + neqns ! create space for diagonals 
      allocate(Ai(nnz))
      allocate(Aj(nnz))
      allocate(Av(nnz))
      nnz = 0
       do i=1, neqns
          do  k=ia(i) , ia(i+1)-1
                 if (i .ne. ja(k)) then
                    nnz = nnz + 1
                    Ai(nnz) = i
                    Aj(nnz) = ja(k)
                    Av(nnz) = -abs(a(k))
                 end if
          end do
      end do
      do i=1, neqns
          do  k=ia(i) , ia(i+1)-1
              if (i .ne. ja(k)) then
                  nnz = nnz + 1
                    Ai(nnz) = ja(k)
                    Aj(nnz) = i
                    Av(nnz) = -abs(a(k))
              end if
          end do
      end do
!     deallocate(B) 
!     deallocate(iB) 
!     deallocate(jB) 
      allocate(D(neqns))
      D = 0.0d0
      do i=1, nnz  
          D(Ai(i)) = D(Ai(i))  +abs(Av(i)) 
      end do
      !nrm = 2*maxval(D)
      f(1,1) = 2*maxval(D) 
      write(*,*) 'min/max on the diagonal :', minval(abs(D)) , maxval(abs(D))  
      !do i=1,neqns 
      !    if (abs(D(i)) .eq. 0.0d0 ) write(*,*) 'row : ' , i , ' is zero.' 
      !end do 
      !!!      perturb = 1.0d-12*f(1,1) 
      perturb = 1.0d-10 * f(1,1) 
      do i =1 , neqns 
         nnz = nnz + 1 
         Ai(nnz) = i 
         Aj(nnz) = i 
         Av(nnz) = D(i)  
      end do 
      !deallocate(D) 

     
    
      allocate(BB(nnz))
      allocate(jBB(nnz))
      allocate(iBB(neqns+1))
      call coocsr(neqns,nnz, Av, Ai, Aj, BB, jBB , iBB)


      deallocate(Av)
      deallocate(Ai)
      deallocate(Aj)
!      nzmax = nnz

           

      
end if 




		t1 = MPI_WTIME()
	
                n = neqns 
		call MPI_BCAST(kk , 1, MPI_INTEGER,0,MPI_COMM_WORLD,error)
		call MPI_BCAST(n , 1, MPI_INTEGER,0,MPI_COMM_WORLD,error)
		call MPI_BCAST(tol , 1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,error)





                if (.not. allocated (displs)) then
	 	 allocate(displs(nb_procs), stat = alloc_error)
	    	 if(alloc_error /= 0) write(*,*)'------------------------ALLOCATION FAILED displs-------------------------------'
                end if 
		

                if (.not. allocated (rcounts)) then
		 allocate(rcounts(nb_procs), stat = alloc_error)
 		 if(alloc_error /= 0) write(*,*)'------------------------ALLOCATION FAILED rcounts-------------------------------'
                end if 





		part1 = n/nb_procs
		part1_LAST = n - (nb_procs -1 )*part1
		do i =1, nb_procs-1
			rcounts(i) = part1
			displs(i) = (i-1)*part1
		end do
		rcounts(nb_procs) = part1_LAST
		displs(nb_procs) = (nb_procs-1)*part1 
                nr = rcounts(rank+1) 
                nloc = nr 

                ratio = min(part1,part1_LAST) 


                if (.not. allocated (p_hat)) then 
                 allocate(p_hat(n,nrhs),stat = alloc_error)
                 p_hat =0.0d0 
                 if(alloc_error /= 0) write(*,*)'------------------------ALLOCATION FAILED p_hat -----------------------------------'
               end if 

               if (.not. allocated (temp)) then
                allocate(temp(n,nrhs),stat = alloc_error) 
                if(alloc_error /= 0) write(*,*)'------------------------ALLOCATION FAILED temp ----------------------------------'
                end if 


t0 = MPI_WTIME()

if (rank .eq. 0 ) then 
       kmax = 0 
       do i =1, n 
          do  k=iBB(i) , iBB(i+1)-1 
             kmax = max(kmax, (jBB(k)-i))      
          end do 
       end do 

!       if (ratio .lt. kmax ) then 
!               write(*,*) '****************************************************************'
!               write(*,*) '*WARNING USING  SLOWER MATVEC: scalability might get affected  *' 
!               write(*,*) '*SOLUTION: use 1 MPI process per node and max number of threads*'
!               write(*,*) '*within the node to enable MPI/OpenMP  parallism  and solve a  *'
!               write(*,*) '*larger system in order to use more nodes                      *'
!               write(*,*) '****************************************************************' 
!               write(*,*) 'kmax:' , kmax, 'nr:', part1,'number of MPI processes:', nb_procs
!               write(*,*) '****************  END OF WARNING MESSAGE  **********************'
!       end if 

end if 



call MPI_BCAST(ratio , 1, MPI_INTEGER,0,MPI_COMM_WORLD,error)
call MPI_BCAST(kmax , 1, MPI_INTEGER,0,MPI_COMM_WORLD,error)

                if (.not. allocated (SEND1)) then
                 allocate(SEND1(kmax,nrhs), stat = alloc_error)
                 if(alloc_error /= 0) write(*,*)'------------------------ALLOCATION FAILED SEND1-------------------------------'
                end if



                if (.not. allocated (SEND2)) then
                 allocate(SEND2(kmax,nrhs), stat = alloc_error)
                 if(alloc_error /= 0) write(*,*)'------------------------ALLOCATION FAILED SEND2-------------------------------'
                end if


                if (.not. allocated (RECV1)) then
                 allocate(RECV1(kmax,nrhs), stat = alloc_error)
                 if(alloc_error /= 0) write(*,*)'------------------------ALLOCATION FAILED RECV1-------------------------------'
                end if


                if (.not. allocated (RECV2)) then
                 allocate(RECV2(kmax,nrhs), stat = alloc_error)
                 if(alloc_error /= 0) write(*,*)'------------------------ALLOCATION FAILED RECV2-------------------------------'
                end if


       !! ALL ranks do together 
       do i =2 , nb_procs
            if (rank .eq. 0 ) nz = iBB(displs(i)+rcounts(i)+1) -iBB(displs(i)+1)
            if (rank .eq. 0) then
                 call MPI_SEND(nz,1,MPI_INTEGER,i-1,100+i,MPI_COMM_WORLD,  error)
            end if
            if (rank .eq. i-1) then
                call MPI_RECV(nz,1,MPI_INTEGER,0,100+i,MPI_COMM_WORLD,status,error)
            end if
       end do
       if (rank .eq. 0) nz = iBB(displs(1)+rcounts(1)+1) - iBB(displs(1)+1)
       nr = rcounts(rank+1)

       if (.not. allocated (Aloc)) then
        allocate(Aloc(nz),stat = alloc_error)
        if(alloc_error /= 0) write(*,*)'------------------------ALLOCATION FAILED Aloc --------------------------------'
       end if 

       if (.not. allocated (jAloc)) then 
        allocate(jAloc(nz),stat = alloc_error)
        if(alloc_error /= 0) write(*,*)'------------------------ALLOCATION FAILED jAloc --------------------------------'
       end if 
 
      
       if (.not. allocated (iAloc)) then      
        allocate(iAloc(nr+1),stat = alloc_error)
        if(alloc_error /= 0) write(*,*)'------------------------ALLOCATION FAILED iAloc --------------------------------' 
       end if 


       do i =2 , nb_procs
           if (rank .eq. 0) then
               nz = iBB(displs(i)+rcounts(i)+1) - iBB(displs(i)+1)
               call MPI_SEND(iBB(displs(i)+1),rcounts(i)+1,MPI_INTEGER,i-1,500+i,MPI_COMM_WORLD,error)
               call MPI_SEND(jBB(iBB(displs(i)+1)),nz,MPI_INTEGER,i-1,600+i,MPI_COMM_WORLD, error)
               call MPI_SEND(BB(iBB(displs(i)+1)),nz,MPI_DOUBLE_PRECISION,i-1,700+i,MPI_COMM_WORLD, error)
           end if
           if (rank .eq. i-1) then
               call MPI_RECV(iAloc(1),rcounts(i)+1,MPI_INTEGER,0,500+i,MPI_COMM_WORLD,status,error)
               idum = iAloc(1)
               do j=1,nr+1
                iAloc(j)  = iAloc(j) - idum  +1
               end do
               call MPI_RECV(jAloc(1),nz,MPI_INTEGER,0,600+i,MPI_COMM_WORLD,status,error)
               call MPI_RECV(Aloc(1),nz,MPI_DOUBLE_PRECISION,0,700+i,MPI_COMM_WORLD,status, error)
            end if
      end do
            if (rank .eq. 0) nz = iBB(displs(1)+rcounts(1)+1) -iBB(displs(1)+1)
            if (rank .eq. 0) then
                do i=1,nr+1
                  iAloc(i) = iBB(i)
                end do
                do i=1,nz
                   jAloc(i) = jBB(i)
                   Aloc(i) = BB(i)
                end do
            end if

      if ( rank .ne. nb_procs-1) then
          i1 =  1 + rank*part1
          i2 =   (rank+1)*part1
       else
          i1 =  1 + rank*part1
          i2 = n
        end if
       

       allocate(Di(nr)) 
!       if (rank.eq.0) write(87,*) D 
       call MPI_SCATTERV(D, rcounts , displs, MPI_DOUBLE_PRECISION,Di,nr , MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,error)
       if (rank .eq. 0 ) deallocate(D) 



!!!!!!!!!! NEW MATVEC 

!!!!!!!!!! NEW MATVEC

     allocate(reqsend(nb_procs))
     allocate(reqrecv(nb_procs))

     allocate(map(nb_procs, n ))
     map = .false.
     do i = 1, nloc
            do k = iAloc(i), iAloc(i+1)-1
!                write(10+rank,*) i, jAloc(k), Aloc(k)
                map(rank+1,jAloc(k)) = .true.
            end do
     end do
     map(rank+1,displs(rank+1)+1:displs(rank+1)+rcounts(rank+1))= .false.
     do i =1 , nb_procs
        call MPI_BCAST(map(i,:) , n, MPI_LOGICAL,i-1,MPI_COMM_WORLD,error)
     end do
     allocate(n_send(nb_procs))
     allocate(n_recv(nb_procs))
     n_send = 0
     n_recv = 0
     do i=1,nb_procs
        do j =displs(i)+1, displs(i)+rcounts(i)
              if (map(rank+1,j))  n_recv(i) = n_recv(i) + 1 ;
        end do
        do j =displs(rank+1)+1, displs(rank+1)+rcounts(rank+1)
              if (map(i,j))  n_send(i) = n_send(i) + 1 ;
        end do
     end do
     allocate(SENDBUF(maxval(n_send)*nrhs*nb_procs))
     allocate(RECVBUF(maxval(n_recv)*nrhs*nb_procs))
     allocate(sendmap(maxval(n_send), nb_procs))
     allocate(recvmap(maxval(n_recv), nb_procs))
     sendmap = 0
     recvmap = 0


     do i=1,nb_procs
        cnt = 0
        do j =displs(i)+1, displs(i)+rcounts(i)
            if (map(rank+1,j)) then
               cnt = cnt +1
               recvmap(cnt,i ) = j
            end if
        end do
        cnt =0
        do j =displs(rank+1)+1, displs(rank+1)+rcounts(rank+1)
             if (map(i,j)) then
              cnt = cnt + 1
              sendmap(cnt,i) = j
             end if
        end do
     end do


end if 

if (step .eq. 999) then  
if (rank .eq. 0) then 
  deallocate(BB) 
  deallocate(iBB) 
  deallocate(jBB) 
end if 
deallocate(displs) 
deallocate(rcounts) 
deallocate(p_hat) 
deallocate(temp) 
deallocate(SEND1) 
deallocate(SEND2) 
deallocate(RECV1) 
deallocate(RECV2) 
deallocate(Aloc) 
deallocate(iAloc) 
deallocate(jAloc) 
deallocate(Di) 
deallocate(reqsend) 
deallocate(reqrecv) 
deallocate(map)
deallocate(n_send)
deallocate(n_recv) 
deallocate(SENDBUF) 
deallocate(RECVBUF)
deallocate(sendmap)
deallocate(recvmap) 
end if 

END SUBROUTINE MATVEC
!        end program main 


   subroutine my_getbwd(n,nnz, ja,ia,kl,ku)
      integer,dimension(1:nnz),intent(in) :: ja
      integer,intent(in) :: n, nnz 
      integer,dimension(1:n+1),intent(in) ::ia
!!!!!!!!!!!!
      integer,intent(out) ::  kl,ku
      integer :: dst,i,k
      kl = - n
      ku = - n
      do  i=1,n
         do  k=ia(i),ia(i+1)-1
            dst = i-ja(k)
            kl = max(kl,dst)
            ku = max(ku,-dst)
         end do
      end do
    end subroutine my_getbwd
