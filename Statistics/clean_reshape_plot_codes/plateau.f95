!*********************************************************************
!Written by Rachith Aiyappa at IISc, Bangalore
!Email rachithaiyappa96@gmail.com for any queries

!This program is to read the ouptut files obtained from the simulation
!Specifically, that of  r^2 which is titled as "r_sq_time_00n.txt" 
!Then, perform a sampling average which evolves with time 
!This the output of the program. It is stored as "r_sq_samp_av.txt"
!Output is <r^2>(t)

!The input/change which needs to be made in the program during
!the run/compile is :
!1. Dimension of x and y
!2. The end of i which depends is equal to "n" in "r_sq_time_00n.txt"
!	In other words, max(i) is the number of input files
!3. j corresponds to the minimum number of columns deduced from 
!	looking into all txt files

!Be careful about the maximum unit number which can be asssigned to a
!particular txt file. Fortran 95 supports unit numbers from 10-99
!*********************************************************************
program plateau

implicit none
character(100) :: filename
integer :: i,j,root_rsq_N_time_file_status
integer :: root_rsq_N_samp_av_file_status
integer :: P2_samp_av_file_status,P2_time_file_status
integer :: av_1_status,av_2_status,av_3_status,av_4_status,av_5_status
integer :: x_status,y_status,z_status
integer :: count,ch
integer :: max_j=1, max_i=100
real :: summ,av_1,av_2,av_3,av_4,av_5

real,dimension(:,:),allocatable :: z,y,x
integer, dimension(100) :: to_change

!opening file to store the time dependent sampling average 
!call system('rm root_rsq_N_samp_av.txt')
!open (unit=10, file="root_rsq_N_samp_av.txt",			&
!iostat=root_rsq_N_samp_av_file_status)
!write(*,*) "root_rsq_N_samp_av_file_status=",			&
!root_rsq_N_samp_av_file_status
open (unit=10, file="P2_samp_av.txt",			&
iostat=P2_samp_av_file_status)
write(*,*) "P2_samp_av_file_status=",			&
P2_samp_av_file_status

!write(*,*) "Enter the number of input files"
!read(*,*) max_i

!write(*,*) "Enter the number of columns to be included from each file"
!read(*,*) max_j

allocate(x(max_i,max_j),stat=x_status)
write(*,*) "x_satus=" ,x_status
allocate(y(max_i,max_j),stat=y_status)
write(*,*) "y_satus=" ,y_status
allocate(z(max_i,max_j),stat=z_status)
write(*,*) "z_satus=" ,z_status

do i = 1,max_i,1
!	write(*,*) i
!	write(filename,'("r_sq_time",i3.3,".txt")') i
!	open (unit=i+10, file=filename,				&
!iostat=root_rsq_N_time_file_status)
!	write(*,*) "root_rsq_N_time_file_status",		&
!root_rsq_N_time_file_status
!end do
	write(filename,'("P",i3.3,".txt")') i
	open (unit=i+10, file=filename,				&
iostat=P2_time_file_status)
	write(*,*) "P2_time_file_status",		&
P2_time_file_status
end do

!x(i,j) : i is the file number, j is the column of the file
!so constant i varying j indicates all the data in a single file
count = 0
do i=11,max_i+10
	!write(*,*) i
	read(i,*) (x(i-10,j),y(i-10,j),z(i-10,j),j=1,max_j)
	j = 1
	do j = 1,max_j	!checking for NaN
		if (z(i-10,j).ne.z(i-10,j)) then
			count = count + 1
			to_change(count) = i-10
			exit
		end if 
	end do
	!read(i+1,*) x,y
	!write(*,*) x,y
	!summ = summ + y
end do
write (*,*) count,to_change

do i = 1,max_i,1
	close(i+10)
end do

do i = 1,count	!rewriting the NaN file
	ch = to_change(i)
!	write(filename,'("r_sq_time",i3.3,".txt")') ch
!	open (unit=ch+10,file=filename,status="replace", 	&
!		action="write",position="rewind",		&
!		iostat=root_rsq_N_time_file_status)
!	write(*,*) 						&
!"root_rsq_N_samp_av_file_status=", root_rsq_N_samp_av_file_status

	write(filename,'("P",i3.3,".txt")') ch
	open (unit=ch+10,file=filename,status="replace", 	&
		action="write",position="rewind",		&
		iostat=P2_time_file_status)
	write(*,*) "P2_time_file_status",P2_time_file_status
	do j = 1,max_j
		write(ch+10,*) x(ch-1,j),y(ch-1,j),z(ch-1,j)
	end do
	close(ch+10)
end do

do i = 1,max_i,1 !reopening all files
	write(*,*) i
!	write(filename,'("r_sq_time",i3.3,".txt")') i
!	open (unit=i+10, file=filename,				&
!		iostat=root_rsq_N_time_file_status)
!	write(*,*) 						&
!"root_rsq_N_samp_av_file_status=", root_rsq_N_samp_av_file_status


	write(filename,'("P",i3.3,".txt")') i
	open (unit=i+10, file=filename,iostat=P2_time_file_status)
	write(*,*) "P2_time_file_status",P2_time_file_status
end do

do i=11,max_i+10 !re-reading new files
	write(*,*) i
	read(i,*) (x(i-10,j),y(i-10,j),z(i-10,j),j=1,max_j)
	!write (*,*) j
	!read(i+1,*) x,y
	!write(*,*) x,y
	!summ = summ + y
end do
write (*,*) count,to_change

do i = 1,max_i,1
	close(i+10)
end do


open (unit=11, file="av_1.txt",iostat=av_1_status)
write(*,*) "av_1_status=",av_1_status
open (unit=12, file="av_2.txt",iostat=av_2_status)
write(*,*) "av_2_status=",av_2_status
open (unit=13, file="av_3.txt",iostat=av_3_status)
write(*,*) "av_3_status=",av_3_status
open (unit=14, file="av_4.txt",iostat=av_4_status)
write(*,*) "av_4_status=",av_4_status
open (unit=15, file="av_5.txt",iostat=av_5_status)
write(*,*) "av_5_status=",av_5_status



!This is to store <r^2>(t) in the file
do j=1,max_j
	av_1 = sum(z(1:(max_i/5),j))/max_i
	av_2 = sum(z(21:(2*(max_i/5)),j))/max_i
	av_3 = sum(z(41:(3*(max_i/5)),j))/max_i
	av_4 = sum(z(61:(4*(max_i/5)),j))/max_i
	av_5 = sum(z(81:(5*(max_i/5)),j))/max_i
	!write(10,*) x(1,j),y(1,j),(av_1+av_2+av_3+av_4+av_5)
	write(10,*) x(1,j),y(1,j),sum(z(1:max_i,j))/max_i
	write(11,*) x(1,j),y(1,j),av_1
	write(12,*) x(1,j),y(1,j),av_2
	write(13,*) x(1,j),y(1,j),av_3
	write(14,*) x(1,j),y(1,j),av_4
	write(15,*) y(1,j),y(1,j),av_5
end do

end program 
