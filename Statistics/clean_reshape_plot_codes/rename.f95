!*********************************************************************
!Written by Rachith Aiyappa at IISc, Bangalore
!Email rachithaiyappa96@gmail.com for any queries

!This script contains code to rename any folder once data is generated
!However, this can also be used to create once folder and move the 
!contents of the other folders into this.
!This was mainly written inorder to copy/move an executable program
!into a directory where the data files exists.
!*********************************************************************


program rename 

implicit none 
character(100) :: folder
character(400) :: re_name
character(200) :: path

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

do f=1,9
	do e=1,9
		write(*,*) "gamma_rand=",gamma_rand(f)
		write(*,*) "gamma_att=",gamma_att(e)
		write(folder,'(a,f4.2,a,f4.2)') 'rand_',	  &
			gamma_rand(f),'_att_',gamma_att(e)
!		write(re_name,'(a,a,a)') 'mv ',		 	  &
!'100_\<r_2\>\(t\).png ', folder
		write(path,'(a)') '/home/arun/rachith/v_3/plots'
		write(re_name,'(a,a,a)') 'cp ',folder, path
		call chdir(folder)
		call system(re_name)
		call chdir('..')
	end do
end do
!******************************************************************
end program rename
