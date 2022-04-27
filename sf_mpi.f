        program Sq_calculator

        include 'mpif.h'
        include 'inc.parameters'
        
C parameter
C a             xyz array of all molecules
C mul           result of q_vec (1 by 3) x a(:,4:)' (3 by num_atoms)

        parameter (c_pi=3.141592)

        integer i, q_1, q_2, q_3
        real q_vec(3)
        real q_and_val(n*n*n,2)
        real q_and_val_tmp(n*n*n,2)

        real a(num_atoms,3)
        real mul(num_atoms), res

C x1 - x3       dummy for reading in xyz
        real x1, x2, x3

C for mpi
        integer send_data_tag, return_data_tag
        parameter (send_data_tag=9005, return_data_tag=9006)
        integer ierr
        integer status(MPI_STATUS_SIZE)
        integer my_id
        integer an_id
        integer num_procs
        integer root_process
        integer start_row
        integer end_row
        integer avg_rows_per_process
        integer num_rows

        call MPI_Init ( ierr )

        open(1,file='input.xyz') ! read in a dump file
        do i = 1, num_atoms, 1
           read(1,*) x1, x2, x3
           a(i,1)=x1
           a(i,2)=x2
           a(i,3)=x3
        enddo 
        close(1)

        root_process=0
        
        call sortrow(a,num_atoms) ! sort rows with the first column ascending

        call MPI_Comm_Rank(MPI_Comm_World,my_id,ierr)
        call MPI_Comm_Size(MPI_Comm_World,num_procs,ierr)

        avg_rows_per_process=(n*n*n)/num_procs

        do an_id = 0, num_procs-1
           if (my_id .eq. an_id) then
              start_row=(an_id*avg_rows_per_process)+1
              end_row=start_row+avg_rows_per_process-1
              if (an_id .eq. (num_procs-1)) end_row=
     &        n*n*n

              num_rows=end_row-start_row+1
           endif
        enddo

        do i=start_row, end_row
           ! vector to each point in the reciprocal space
           q_1=(i-1)/n**2
           q_2=mod(i-1,n**2)/n
           q_3=mod(mod(i-1,n**2),n)
           q_vec=(/q_1, q_2, q_3/)*2*c_pi/L ! delta_q=2*pi/L

           q_modulus=NORM2(q_vec)
           q_and_val_tmp(i,1)=q_modulus

           mul=matmul(a,q_vec)
           if (i .eq. 1) then
              write(*,*) q_vec
           endif

           res=sum(cos(mul))**2+sum(sin(mul))**2
           if (i .eq. 1) then
              write(*,*) res/num_atoms
           endif

           q_and_val_tmp(i,2)=res/num_atoms

        enddo

        if (my_id .eq. 0) then
           
           q_and_val(start_row:end_row,:)=
     &     q_and_val_tmp(start_row:end_row,:)

           do an_id=1,num_procs-1
              start_row=(an_id*avg_rows_per_process)+1
              end_row=start_row+avg_rows_per_process-1
              if (an_id .eq. (num_procs-1)) end_row=
     &        n*n*n

              num_rows=end_row-start_row+1

              call MPI_Recv(q_and_val_tmp(1:num_rows,:),
     &        num_rows*2,MPI_Real,an_id,
     &        MPI_Any_Tag,MPI_Comm_World,status,ierr)

              q_and_val(start_row:end_row,:)=
     &        q_and_val_tmp(1:num_rows,:)
           enddo
        else
           call MPI_Send(q_and_val_tmp(start_row:end_row,:),
     &     num_rows*2,MPI_Real,root_process,
     &     return_data_tag,MPI_Comm_World,ierr)
        endif

        call MPI_Finalize(ierr)

        if (my_id .eq. 0) then
           open(2, file='sq.txt')
           do i=1,n**3
              write (2,*) q_and_val(i,:)
           enddo
           close(2)
        endif
       
        stop
        end

C -----------------------------------------------------------

        subroutine sortrow(A,num_rows)
        real A(num_rows,6), buf(6)
        integer irow, krow

        do irow = 1, num_rows, 1
           krow=minloc(A(irow:num_rows,1),dim=1)+irow-1
           buf(:)=A(irow,:)
           A(irow,:)=A(krow,:)
           A(krow,:)=buf(:)
        enddo

        return
        end

C -----------------------------------------------------------
