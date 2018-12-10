!*********************************************************************
!Written by Rachith Aiyappa at IISc, Bangalore
!Email rachithaiyappa96@gmail.com for any queries

!This enters each of the sudirectories containing data, then collects
!the already existing average and puts all such averages from different
!data subdirectories into one file in the main directory.

!This one file can then be later used for plotting using gnuplot
!*********************************************************************

program onefile
implicit none 

character(100) ::folder
integer :: i,e,f
integer :: root_rsq_N_file_status,root_rsq_N_samp_av_file_status
integer :: P2_file_status,P2_samp_av_file_status
real :: x,y,z

real,dimension(30) :: gamma_att 
real,dimension(30) :: gamma_rand

!opening file to store all P2
!open (unit=10, file="root_rsq_N.txt",iostat=root_rsq_N_file_status)
!write(*,*) "root_rsq_N_file_status=",root_rsq_N_file_status
open (unit=10, file="P2.txt",iostat=P2_file_status)
write(*,*) "P2_file_status=",P2_file_status

!initialising gamma's-------------------------------------------------
gamma_att(1) = 0.15
!write(*,*) gamma_att(1)
gamma_rand(1) = 0.15
do i=2,8,1
	gamma_att(i) = gamma_att(i-1)+0.05
	write(*,*)gamma_att(i)
end do
do i=2,6,1
	gamma_rand(i) = gamma_rand(i-1)+0.05
	write(*,*)gamma_rand(i)
end do
!******************************************************************

do f=1,6
	do e=1,8
!		write(folder,'(a,f4.2,a,f4.2,a)') 'rand_',	   &
!			gamma_rand(f),'_att_',gamma_att(e),        &
!			'/sqrt_r_sq_by_N'
!		write(*,*) folder
!		call chdir(folder)
!		open(unit=11,file="root_rsq_N_samp_av.txt",	   &
!		iostat=root_rsq_N_samp_av_file_status)
!		write(*,*) root_rsq_N_samp_av_file_status

		write(folder,'(a,f4.2,a,f4.2,a)') 'rand_',	   &
			gamma_rand(f),'_att_',gamma_att(e),        &
			'/P2'
		write(*,*) folder
		call chdir(folder)
		open(unit=11,file="P2_samp_av.txt",	   &
		iostat=P2_samp_av_file_status)
		write(*,*) P2_samp_av_file_status

		read(11,*) x,y,z
		write(*,*) x,y,z
		write(10,*) x,y,z
		call chdir('../..')
	end do
end do
!******************************************************************
end program onefile
