program kursach
 include 'mpif.h'
  real, allocatable :: A(:,:), B(:,:), C(:), multimatrix(:,:),partA(:,:), endA(:,:), recvA(:,:), testA(:,:),coord(:)
  integer :: i,j,razB,mperr, razmer,k,l, n, rank, numberline, z,ll, stat(MPI_STATUS_SIZE)
  integer, allocatable :: displ(:), numbersend(:)
  double precision :: Time 
  real :: sumnum=0, start_time, finish_time
  open(1, file="reult.txt")
  n=10000
  razB=2
  !create a matrix  A and B
  allocate (A(n,n),B(razB,razB),endA(n*razB,n*razB), recvA(n*razB,n*razB),coord(n*razB))
  call init_random_seed(seed)
  call random_number(A)
  call random_number(B)
   
  call MPI_INIT(mperr)  !start parallel part 
  call MPI_COMM_SIZE(MPI_COMM_WORLD, razmer, mperr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mperr)
 ! call cpu_time(start_time)
  if (rank .eq. 0) then
     Time=MPI_WTIME()
  end if

  allocate (numbersend(razmer), displ(razmer))
  call MPI_BCAST(razB, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr)
  call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr)
  call MPI_BCAST(B, razB*razB, MPI_REAL, 0, MPI_COMM_WORLD, mperr)
  do i=1, razmer
    numbersend(i)=n/razmer
  end do
  if (rank .eq. 0) then
  j=mod(n,razmer)

  do i=1, j
    numbersend(i)=numbersend(i)+1
  end do
  displ(1)=0
  do i=2, razmer
    displ(i)=displ(i-1)+numbersend(i-1)
  end do
  displ=displ*n
  !print *,displ
  numbersend=numbersend*n
  
  end if
  call MPI_SCATTER(numbersend/n, 1, MPI_REAL, numberline, 1, MPI_REAL, 0, MPI_COMM_WORLD, mperr) 
  call MPI_BARRIER(MPI_COMM_WORLD)


  allocate(partA(numberline,n))
  


  call MPI_SCATTERV(A, numbersend, displ, MPI_REAL, partA, numberline*n , MPI_REAL, 0, MPI_COMM_WORLD, mperr)

       allocate(multimatrix(numberline*razB,n*razB), testA(numberline,n*razB))
   
   do k=1, razB
      do i=1, numberline
         do m=1, razB
            do j=1, n
               multimatrix(i+(k-1)*numberline,j+(m-1)*n)=partA(i,j)*B(k,m)
           end do
       end do
     end do
  end do
 ! do i=1,
  do i=1, razB*numberline
 ! print *,rank, multimatrix(i,:)
  end do
  multimatrix=transpose(multimatrix)
  call MPI_GATHERV(multimatrix, numberline*razB*n*razB, MPI_REAL, endA, numbersend*razB*razB, displ*razB*razB, MPI_REAL, &
  0, MPI_COMM_WORLD, mperr)
  
  if (rank .eq. 0) then 
     endA=transpose(endA)
     !print *,numbersend/n
           !do ll=1,razB*n
     ll=1
     do l=1,razB
        sumnum=0
        do k=1,razmer
           do i=1,numbersend(k)/n
              recvA(ll,:)=endA(i+sumnum+(l-1)*numbersend(k)/n,:)
          !    coord(i+sumnum+(l-1)*n)=(k-1)
              !recvA(i+sumnum+(l-1)*n,:)=endA(ll,:)
              ll=ll+1
           end do
           sumnum=sumnum+numbersend(k)/n*razB 
        end do
     end do
  
   !do i=1, n*razB
   !   write (1, *)endA(i,:)
   !end do
   !   write(1, *)" "
   !do i=1, n*razB
   !   write (1, *)recvA(i,:)
   !end do
   Time=MPI_WTIME()-Time
   print *,Time, rank
   call cpu_time(finish_time)
   print *, 'Время вычислений = ', finish_time - start_time
   end if
   

   
  call MPI_FINALIZE(mperr)

end



subroutine init_random_seed()
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed

  call random_seed(size = n)
  allocate(seed(n))

  call system_clock(count=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)

  deallocate(seed)
end subroutine init_random_seed

