!*********************************************************************
!Written by Rachith Aiyappa at IISc, Bangalore
!Email rachithaiyappa96@gmail.com for any queries

!This script copies an executable to each subdirectory. Then this same 
!script runs the executable in each subdirectory and perfroms an 
!average of the data present accross multiple files present in the 
!sudirectory
!*********************************************************************

program samp
implicit none 
character(100) :: folder
character(200) :: copy

real,dimension(30) :: gamma_att 
real,dimension(30) :: gamma_rand

integer :: i,e,f


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

!copying executable plateau into all-----------------------------------
do f=1,6
	do e=1,8
!		write(folder,'(a,f4.2,a,f4.2,a)') 'rand_',	   &
!			gamma_rand(f),'_att_',gamma_att(e),	   &
!'/sqrt_r_sq_by_N'
write(folder,'(a,f4.2,a,f4.2,a)') 'rand_',	   &
			gamma_rand(f),'_att_',gamma_att(e),	   &
'/P2'
		write(*,*) folder
		write(copy,'(a,a)') 'cp plateau ',folder
		write(*,*) copy
		call system(copy)
	end do
end do
!******************************************************************

!loop to run plateau in all---------------------------------------
do f=1,6
	do e=1,8
!	write(*,*) "gamma_rand=",gamma_rand(f)
!		write(*,*) "gamma_att=",gamma_att(e)
!		write(folder,'(a,f4.2,a,f4.2,a)') 'rand_',	  &
!			gamma_rand(f),'_att_',gamma_att(e),	  &
!'/sqrt_r_sq_by_N'
		write(folder,'(a,f4.2,a,f4.2,a)') 'rand_',	  &
			gamma_rand(f),'_att_',gamma_att(e),	  &
'/P2'
		call chdir(folder)
		call system('./plateau')
		call chdir('../..')
	end do
end do
!******************************************************************

end program samp
