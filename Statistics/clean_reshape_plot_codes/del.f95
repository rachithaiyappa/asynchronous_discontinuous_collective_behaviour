!*********************************************************************
!Written by Rachith Aiyappa at IISc, Bangalore
!Email rachithaiyappa96@gmail.com for any queries

!This script is used to delete any contents of the subdirectory if
!required as per user requirement
!*********************************************************************

program del
implicit none
character(100) :: folder

real,dimension(30) :: gamma_att 
real,dimension(30) :: gamma_rand

integer :: i,e,f


!initialising gamma's--------------------------------------------------
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
!*********************************************************************

!deleting executable plateau in all-----------------------------------
do f=1,6
	do e=1,8
		write(folder,'(a,f4.2,a,f4.2,a)') 'rand_',	  &
			gamma_rand(f),'_att_',gamma_att(e),'/P2'
		write(*,*) folder
		call chdir(folder)
		call system('pwd')
		call system('rm plateau')
		call system('rm r_sq_samp_av.txt')
		call system('rm av_1.txt')
		call system('rm av_2.txt')
		call system('rm av_3.txt')
		call system('rm av_4.txt')
		call system('rm av_5.txt')
		call chdir('../..')
	end do
end do
!*********************************************************************
end program del
