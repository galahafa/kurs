program kursach
   ! This program make Kronecker product with using MPI technology
   include 'mpif.h'
   real, allocatable :: A(:,:), B(:,:), C(:), multimatrix(:,:),partA(:,:), endA(:,:), recvA(:,:)
   !A and B -initial matrix; partA, endA-matrix using for send and recv matrix; multimatrix-result multiply matrix on each thread
   ! recvA - result
   integer :: i,j,razB,mperr, razmer,k,l, n, rank, numberline, z,ll
   !razB - size matrix B, n - size matrix A, rank - personal number each threads, razmer -number of threads
   !numberline - number line which get thread, other varaible using for count
   integer, allocatable :: displ(:), numbersend(:)
   !displ and numbersend using for MPI_SCATTERV 
   double precision :: Time 
   real :: sumnum=0, start_time, finish_time
   open(1, file="result.txt")
   !define size of matrix A and B
   n=4
   razB=2
   !create a matrix  A and B
   allocate (A(n,n),B(razB,razB),endA(n*razB,n*razB), recvA(n*razB,n*razB))
   call init_random_seed(seed)
   call random_number(A)
   call random_number(B)
   ! start main part 
   call MPI_INIT(mperr)  !start parallel part 
   call MPI_COMM_SIZE(MPI_COMM_WORLD, razmer, mperr) !define number of threads
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mperr) !give each thread the number
   call cpu_time(start_time) !start count cpu_time
   if (rank .eq. 0) then    !start count program time
      Time=MPI_WTIME()
   end if
   allocate (numbersend(razmer), displ(razmer))
   !send size of matrix A and B, and matrix B to all threads
   call MPI_BCAST(razB, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr)
   call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr)
   call MPI_BCAST(B, razB*razB, MPI_REAL, 0, MPI_COMM_WORLD, mperr)
   ! prepearing to send matrix A
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
      numbersend=numbersend*n
   end if
   call MPI_SCATTER(numbersend/n, 1, MPI_REAL, numberline, 1, MPI_REAL, 0, MPI_COMM_WORLD, mperr) 
   call MPI_BARRIER(MPI_COMM_WORLD)
   allocate(partA(numberline,n))
   partA=0
   A=transpose(A)
   ! send matrix part of matrix A to all threads
   call MPI_SCATTERV(A, numbersend, displ, MPI_REAL, partA, numberline*n , MPI_REAL, 0, MPI_COMM_WORLD, mperr)
   allocate(multimatrix(numberline*razB,n*razB))
   !multiply part matrix A and B on all threads
   do k=1, razB
      do i=1, numberline
         do m=1, razB
            do j=1, n
               multimatrix(i+(k-1)*numberline,j+(m-1)*n)=partA(i,j)*B(k,m)
            end do
         end do
      end do
   end do
   !collect result of multiply on each process 
   multimatrix=transpose(multimatrix)
   call MPI_GATHERV(multimatrix, numberline*razB*n*razB, MPI_REAL, endA, numbersend*razB*razB, displ*razB*razB, MPI_REAL, &
   0, MPI_COMM_WORLD, mperr)
  
   if (rank .eq. 0) then 
      endA=transpose(endA)
      ll=1
      !change incorrect sequence of line after recieve
      do l=1,razB
         sumnum=0
         do k=1,razmer
            do i=1,numbersend(k)/n
               recvA(ll,:)=endA(i+sumnum+(l-1)*numbersend(k)/n,:)
               ll=ll+1
            end do
            sumnum=sumnum+numbersend(k)/n*razB 
         end do
      end do
      !uncomment this part if you nead write result to file
      !do i=1, n*razB
      !  write (1, *)recvA(i,:)
      !end do
      Time=MPI_WTIME()-Time
      print *,Time, rank
      call cpu_time(finish_time)
      print *, 'Время вычислений = ', finish_time - start_time
   end if
   call MPI_FINALIZE(mperr)

end





