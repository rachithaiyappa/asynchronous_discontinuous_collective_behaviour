!*********************************************************************
!Written by Rachith Aiyappa at IISc, Bangalore
!Email rachithaiyappa96@gmail.com for any queries

!This program plots <r^2(t)> for each value of gamma_rand and gamm_att
!Additionally it also copies command.plt into all folders
!*********************************************************************
program plot_trans
implicit none 
character(100) :: folder
character(200) :: copy

real,dimension(30) :: gamma_att 
real,dimension(30) :: gamma_rand

integer ::  i,e,f

!initialising gamma's-------------------------------------------------
gamma_att(1) = 0.1
!write(*,*) gamma_att(1)
gamma_rand(1) = 0
do i=2,9,1
	gamma_att(i) = gamma_att(i-1)+0.05
	write(*,*)gamma_att(i)
end do
do i=2,9,1
	gamma_rand(i) = gamma_rand(i-1)+0.05
	write(*,*)gamma_rand(i)
end do
!******************************************************************

!copying gnu coomands into all-----------------------------------
do f=1,9
	do e=1,9
		write(folder,'(a,f4.2,a,f4.2)') 'rand_',	  &
			gamma_rand(f),'_att_',gamma_att(e)
		write(copy,'(a,a)') 'cp command.plt ',folder !copy gnu plot commands
		call system(copy)
	end do
end do
!*****************************************************************

!loop to run gnuplot in all-------------------------------------------
do f=1,9
	do e=1,9
		write(*,*) "gamma_rand=",gamma_rand(f)
		write(*,*) "gamma_att=",gamma_att(e)
		write(folder,'(a,f4.2,a,f4.2)') 'rand_',	  &
			gamma_rand(f),'_att_',gamma_att(e)      
		call chdir(folder)
		call system ('gnuplot -p command.plt')
		call chdir('..')
	end do
end do
!******************************************************************

end program plot_trans
