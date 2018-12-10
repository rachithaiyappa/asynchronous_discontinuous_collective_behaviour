!*********************************************************************
!Written by Rachith Aiyappa at IISc, Bangalore
!Email rachithaiyappa96@gmail.com for any queries
!This contains the main code which models the collective behaviour
!of N fish without a boundary
!v_3_2 is aimed at seeing the time evolution of radius of gyration
!Like v_2_1, this does not contain information at equal time intervals
!Hence interpolation is required to achive this

!r_sq(t) = (sum_{i,j}(r_i(t)-r_j(t))^2)/(0.5*N(N-1)) 

!(For personal reference) : Another change being the intial postions 
!of the fish is now made
!random

!This contains the loop for running for multiple values of gamma_att 
!and gamma_rand. 
!Minor changes can be implemented to the code inorder to make gamma_ali
!a variable too. At the moment, gamma_ali is set to zero.

!The data of the radius of gyration is required to be plotted 
!against time in order to identify the bounded group size and hence
!the transient time
!*********************************************************************

program nfish

use :: fish_func

implicit none
!Declaration of variable----------------------------------------------

character(100) :: filename
character(100) :: cmd

real, parameter :: pi=3.14159265
!real, parameter :: gamma_att = 0.55
real, parameter :: gamma_ali = 0
!real, parameter :: gamma_rand = 0.04
integer, parameter :: max_size=1000000
real, parameter :: tau_0 = 1.6

integer :: time_limit
!integer :: size_xy=10
integer :: N, i, j, a, c
integer :: e,f
integer :: x_ini_status,y_ini_status, phi_ini_status
integer :: x_status,y_status,phi_status
integer :: x_temp_status,y_temp_status
integer :: d_status,v_status,tau_status,l_status
integer :: time_status,time_temp_status,time_mod_status,cont_t_status
integer :: phi_focal_status,phi_other_status,delta_phi_status
integer :: traj_file_status,fnn_file_status,pol_file_status
integer :: r_sq_file_status
integer :: pol_status,r_sq_status,d_fnn_status
integer :: psi_status
integer :: out_dist_calc_status
integer :: col,count_number
integer :: b_status,noise_status
integer :: chdir_status,md_status
integer :: columns_status,column_file_status
integer :: idx,idx_new,idx_prev
integer :: samp_no

!real :: gamma_att,gamma_ali,gamma_rand
real :: del_phi_att_total,del_phi_ali_total
real :: d_first_nn_all
real :: focal_x, focal_y
real :: x_first_nn,y_first_nn
real :: val
real :: x_comp,y_comp
real :: del_phi_ali_tot, del_phi_att_tot
real :: r
real :: cos_term, sin_term, pol_sum, pol_mod
real :: summ

!Static arrays
integer, dimension(max_size) :: fish_prev,fish_new
integer, dimension(max_size) :: foc_x,foc_y
real,dimension(30) :: gamma_att 
real,dimension(30) :: gamma_rand

!Dynamically allocatable arrays
integer, dimension(:), allocatable :: columns
real, dimension(:),allocatable :: x_ini,y_ini,phi_ini
real, dimension(:),allocatable :: pol,r_sq
real, dimension(:,:),allocatable :: x,y,phi
real, dimension(:,:),allocatable :: x_temp,y_temp
real, dimension(:,:),allocatable :: d,out_dist_calc
real, dimension(:,:),allocatable :: v,tau,l,out_v,out_tau
real, dimension(:,:),allocatable :: time,time_temp,time_mod,cont_t
real, dimension(:,:),allocatable :: b
real, dimension(:,:),allocatable :: phi_other,phi_focal,delta_phi
real, dimension(:,:),allocatable :: psi
real, dimension(:,:),allocatable :: g
real, dimension(:,:),allocatable :: d_fnn

!*********************************************************************
!opening the file to store the trajectories
!open (unit = 10,file="traj.txt",status="replace", 	&
!		action="write",position="rewind",	&
!		iostat=traj_file_status)
!write(*,*) "traj_status=", traj_file_status

!opening the file to store fnn distance
!open (unit = 11, file="dist_fnn.txt",status="replace",	&
!		action="write",position="rewind",	&
!		iostat=fnn_file_status)
!write(*,*) "fnn_file_status",fnn_file_status

!opening the file to store polarisation P2
!open (unit = 12, file="P2.txt",status="replace",	&
!		action="write",position="rewind",	&
!		iostat=pol_file_status)
!write(*,*) "p2_file_status",pol_file_status

!opening the file to store time evo of square of radius of gyration
!write(filename,'("myfile",i3.3,".txt")') i
!open (unit=13, file="r_sq_time.txt",status="replace",	&
!		action="write",position="rewind",	&
!		iostat=r_sq_file_status)
!write(*,*) "r_sq_file_status",r_sq_file_status

!*********************************************************************

!User input-----------------------------------------------------------
write(*,*) "Enter the number of fish"
read(*,*) N
write(*,*) "Enter the time of simulation. Value entered = 2*time(s)"
read(*,*) time_limit
write(*,*) "Enter the number of samples required"
read(*,*) samp_no
!*********************************************************************

!Allocating initial conditions----------------------------------------
allocate(x_ini(N),stat = x_ini_status)
allocate(y_ini(N),stat = y_ini_status)
allocate(phi_ini(N),stat = phi_ini_status)

write(*,*) "x_ini_status=", x_ini_status !check if dyn mem allo occured
write(*,*) "y_ini_status=",y_ini_status	!check if dyn mem allo occured
write(*,*) "phi_ini_status=",phi_ini_status !check if dyn mem allo occured
!*********************************************************************

!Allocating x,y,phi for computation-----------------------------------
allocate(x(N,max_size),stat = x_status)
allocate(y(N,max_size),stat = y_status)
allocate(phi(N,max_size),stat = phi_status)

write(*,*) "x_status=",x_status	!check if dyn mem allo occured
write(*,*) "y_status=",y_status	!check if dyn mem allo occured
write(*,*) "phi_status=",phi_status !check if dyn mem allo occured
!*********************************************************************

!---------------------------------------------------------------------
!Allocating x_temp, y_temp for trajectories and movies----------------
allocate(x_temp(N,max_size),stat = x_temp_status)
allocate(y_temp(N,max_size),stat = y_temp_status)

write(*,*) "x_temp_status=",x_temp_status !check if dyn mem allo occured
write(*,*) "y_temp_status=",y_temp_status !check if dyn mem alloo ccured
!*********************************************************************

!Allocating time------------------------------------------------------
allocate(time(N,max_size),stat=time_status)
allocate(time_temp(N,max_size),stat=time_temp_status)
allocate(time_mod(N,max_size),stat=time_mod_status)
allocate(cont_t(N,max_size),stat=cont_t_status)

write(*,*) "time_status=",time_status
write(*,*) "time_temp_status=",time_temp_status
write(*,*) "time_mod_status=",time_mod_status
write(*,*) "cont_t_status=",cont_t_status
!*********************************************************************

!allocating d,l,v,tau-------------------------------------------------
allocate(d(N,max_size),stat=d_status)
allocate(v(N,max_size),stat=v_status)
allocate(tau(N,max_size),stat=tau_status)
allocate(l(N,max_size),stat=l_status)

write(*,*) "d_status=",d_status	
write(*,*) "v_status=",v_status
write(*,*) "tau_status=",tau_status
write(*,*) "l_status=",l_status	
!*********************************************************************

!allocation of d_fnn--------------------------------------------------
allocate(d_fnn(N,max_size),stat=d_fnn_status)
write(*,*) "d_fnn_status=",d_fnn_status
!*********************************************************************

!delta_phi and corresponding allocations------------------------------
allocate(phi_focal(N,max_size),stat=phi_focal_status)
write(*,*) "phi_focal_status=",phi_focal_status

allocate(phi_other(N,max_size),stat=phi_other_status)
write(*,*) "phi_other_status=",phi_other_status

allocate(delta_phi(N,max_size),stat=delta_phi_status)
write(*,*) "delta_phi_status=",delta_phi_status
!*********************************************************************

!allocation of psi----------------------------------------------------
allocate(psi(N,max_size),stat=psi_status)
write(*,*) "psi_status=", psi_status
!*********************************************************************

!allocation of stats--------------------------------------------------
allocate(pol(max_size),stat = pol_status)
write(*,*) "pol_status=", pol_status

allocate(r_sq(max_size),stat = r_sq_status)
write(*,*) "r_sq_status=", r_sq_status
!*********************************************************************

!allocating to keep track of minimum number of columns in the samples-
allocate(columns(samp_no),stat=columns_status)
write(*,*) "columns_status=",columns_status
!*********************************************************************

!Miscellaneous allocation---------------------------------------------
allocate(b(N,max_size),stat=b_status)
write(*,*) "b_status=",b_status

allocate(g(N,max_size),stat=noise_status)
write(*,*) "noise_status=",noise_status

allocate(out_dist_calc(N,1),stat = out_dist_calc_status)
write(*,*) "out_dist_clac=",out_dist_calc_status
!*********************************************************************

!initialising gamma's--------------------------------------------------
gamma_att(1) = 0.10
!gamma_att(2) = 0.15
!write(*,*) gamma_att(1)
gamma_rand(1) = 0.40
!do i=2,7,1
!	gamma_att(i) = gamma_att(i-1)+0.05
!	write(*,*)gamma_att(i)
!end do
!do i=2,3,1
!	gamma_rand(i) = gamma_rand(i-1)+0.05
!	write(*,*)gamma_rand(i)
!end do
!*********************************************************************


!starting loop for gamma_att and noise--------------------------------
do f=1,1
	do e=1,1
		write(*,*) "gamma_rand=",gamma_rand(f)
		write(*,*) "gamma_att=",gamma_att(e)

!creating a new directory to store samples for each (rand,att)---------
write(cmd,'(a,f4.2,a,f4.2)') 'mkdir rand_',gamma_rand(f),	     &
'_att_',gamma_att(e)
write(*,*) cmd
call system(cmd,status=md_status)
!call system("mkdir direc")

write(*,*) "md_status=",md_status
write(cmd,'(a,f4.2,a,f4.2)') 'rand_',gamma_rand(f),	     	     &
'_att_',gamma_att(e)
!call chdir(cmd,status=chdir_status)
call chdir(cmd)
!write(*,*) "chdir_status=",chdir_status


!opening the file to store r^2(t)-------------------------------------
do i = 1,samp_no,1
write(*,*) i
write(filename,'("r_sq_time",i3.3,".txt")') i
open (unit=13, file=filename,status="replace",			     &
		action="write",position="rewind",		     &
		iostat=r_sq_file_status)
write(*,*) "r_sq_file_status",r_sq_file_status
!*********************************************************************

!*********************************************************************
call init_random_seed
!*********************************************************************

!initial postions and angles of the fish------------------------------
do j =1,N,1
	call random_number(r)
	x_ini(j) = 10*r
	y_ini(j) = 10*r
	phi_ini(j) = 2.0*pi*r
end do
!write(*,*) x_ini,y_ini,phi_ini	!chk if ini pos of the fish are stored
!*********************************************************************

!copying intial postions to x,y,phi-----------------------------------
do j =1,N,1
	x(j,1) = x_ini(j)
	!write(*,*) x
	y(j,1) = y_ini(j)
	phi(j,1) = phi_ini(j)
end do
!*********************************************************************

!x_temp and y_temp stores info for traj plots-------------------------
x_temp(1:N,1) = x(:,1)
!write(*,*) x_temp(1:N,1)
y_temp(1:N,1) = y(:,1)
!*********************************************************************

!absolute time--------------------------------------------------------
time(:,1) = 0.0
!write(*,*) time(:,1)
!*********************************************************************

!calculating initial distance from focal fish-------------------------
focal_x = minval(x(:,1))
focal_y = minval(y(:,1))
col = 1	

out_dist_calc = dist_calc(col,focal_x,focal_y,x(:,col),y(:,col),N)
d(:,col) = out_dist_calc(:,col)
!write(*,*) d(:,col)
!*********************************************************************

!generating first column of noise-------------------------------------
do a=1,N,1
	g(a,1) = gaus_rand()
end do
!write(*,*) g(:,1)
!*********************************************************************

!generating the first column of velocities----------------------------
col = 1
a = 1
do a=1,N,1
	call random_number(r)
	v(a,col) = generate_v(r)
end do
!write (*,*) v(:,1)
!out_v = generate_v(col,N)
!v(:,col) = out_v(:,col)
!write (*,*) v(:,1)
!out_v = generate_v(col+1,N)
!v(:,col+1) = out_v(:,col+1)
!write(*,*) v(:,1), v(:,2), v(:,3)
!*********************************************************************

!generating the first column of kicking times-------------------------
col = 1
a = 1
do a=1,N,1
	call random_number(r)
	!write (*,*) r
	tau(a,col) = generate_tau(r)
end do
!write (*,*) tau(:,1)
!out_tau = generate_tau(col,N)
!tau(:,col) = out_tau(:,col)
!write (*,*) tau(:,1)
!*********************************************************************

!The next step at which parameters need to be updated
!once the focal fish has been identified
time(:,2) = time(:,1) + tau(:,1)
!write (*,*) time(:,1),time(:,2)
!*********************************************************************

!Storing the first column of time temp--------------------------------
idx = minloc(time(:,2),dim=1)
time_temp(idx,1) = time(idx,2)
write (*,*) "blah"
!write(*,*) idx,time_temp(:,1)
a = 1
do while (a .le. N)
	if (a .ne. idx) then
		time_temp(a,1) = time(a,1)
	end if
	a = a + 1
end do
!write (*,*) idx,time_temp(:,1)
!*********************************************************************

summ = 0.0
do a=1,N,1
	do c=a+1,N,1
		summ= summ+(((x_temp(a,col)-x_temp(c,col))**2)	     &
		+((y_temp(a,col)-y_temp(c,col))**2))
	end do
end do
	r_sq(col) = (2.0/(N*(N-1)))*summ
	write(13,*) time_mod(1,col),r_sq(col)

!Variable defined inorder to enter the loop---------------------------
count_number = 0
col = 2
!*********************************************************************

do while (time_limit > maxval(time(:,col)))
	!write(*,*) "max_time=",maxval(time(:,col))
	count_number = count_number+1

	!index of the previous focal fish	
	idx_prev = minloc(time(:,col-1),dim=1)
	fish_prev(col) = idx_prev

	!index of the new focal fish
	idx_new = minloc(time(:,col),dim=1)
	fish_new(col) = idx_new
	!write(*,*) idx_new

	!storing the same time in all rows of a particular column in
	!order to bring all fish to that time
	time_mod(:,col)=time(idx_new,col)

	!updating time_temp which is used to get cont_t for other fish
	time_temp(idx_new,col)=minval(time(:,col))
	a = 1
	do while (a .le. N)
		if (a .ne. idx_new) then
			time_temp(a,col) = time_temp(a,col-1)
		end if
		a = a + 1
	end do

	!length tranversed by ff before choosing new v,tau
	l(idx_new,col-1) = v(idx_new,col-1)*tau_0*		     &
			(1-(exp(-(tau(idx_new,col-1)/tau_0))))
	
	!bringing ff till it meets its own tau
	cont_t(idx_new,col-1) = tau(idx_new,col-1)
	x(idx_new,col) = x(idx_new,col-1) + ((l(idx_new,col-1)*      &
			((1-exp(-cont_t(idx_new,col-1)/tau_0))/	     &
			(1-exp(-tau(idx_new,col-1)/tau_0)))*         &
			cos(phi(idx_new,col-1))))
	y(idx_new,col) = y(idx_new,col-1) + ((l(idx_new,col-1)*      &
			((1-exp(-cont_t(idx_new,col-1)/tau_0))/	     &
			(1-exp(-tau(idx_new,col-1)/tau_0)))*         &
			sin(phi(idx_new,col-1))))
	b(idx_new,col-1) = ((1-exp(-cont_t(idx_new,col-1)/tau_0))/   &
			(1-exp(-tau(idx_new,col-1)/tau_0)))
	!write(*,*) b(idx_new,col-1)
	foc_x(col) = x(idx_new,col)
	foc_y(col) = y(idx_new,col)

	!length traveresed by other fish to reach pt to give info to ff
	a = 1
	do while (a .le. N)
		if (a .ne. idx_new) then
			l(a,col-1) = v(a,col-1)*tau_0*         	     &
				(1-(exp(-(tau(a,col-1)/tau_0))))
		end if 
		a = a + 1
	end do


	!bringing other fish to meet ff tau
	a = 1
	do while (a .le. N)
		if (a .ne. idx_new) then
			cont_t(a,col-1) = time_temp(idx_new,col) -   &
						time_temp(a,col)
			x(a,col) = x(a,col-1) + ((l(a,col-1)*        &
				((1-exp(-cont_t(a,col-1)/tau_0))/    &
				(1-exp(-tau(a,col-1)/tau_0)))*	     &
				cos(phi(a,col-1))))
			y(a,col) = y(a,col-1) + ((l(a,col-1)*	     &
				((1-exp(-cont_t(a,col-1)/tau_0))/    &
				(1-exp(-tau(a,col-1)/tau_0)))*	     &
				sin(phi(a,col-1))))
			b(a,col-1)=((1-exp(-cont_t(a,col-1)/tau_0))/ &
				   (1-exp(-tau(a,col-1)/tau_0)))       
		end if
		a = a + 1
	end do

	!this has info abt pos of other fish at the kick time of ff
	x_temp(:,col) = x(:,col)
	y_temp(:,col) = y(:,col)
	!write(*,*) y_temp(:,col-1)

	!this has info of updates of only the ff
	!the other fish remain where they were
	a = 1
	do while (a .le. N)
		if (a .ne. idx_new) then
			x(a,col) = x(a,col-1)
			y(a,col) = y(a,col-1)
		end if 
		a = a + 1
	end do

	!calculating the distance from ff to other fish
	!taking into account the ff and its first nearest neighbour(nn)
	focal_x = x_temp(idx_new,col)
	focal_y = y_temp(idx_new,col)
	!write(*,*) x_temp(1:N,col)
	out_dist_calc = dist_calc(col,focal_x,focal_y,		     &
		x_temp(1:N,col),y_temp(1:N,col),N)
	d(:,col) = out_dist_calc(:,1)
	!write(*,*) d(1:N,col)
	
	!calculating fnn distance from all fish for stats
	a = 1
	do while (a .le. N)
		out_dist_calc = dist_calc(col,x_temp(a,col),	     &
			y_temp(a,col),x_temp(:,col),y_temp(:,col),N)
		!dist of other fish from fish a at current time
		idx = minloc(out_dist_calc(:,1),dim=1,mask=	     &
		out_dist_calc(:,1) .gt. minval(out_dist_calc(:,1)))
		x_first_nn = x_temp(idx,col)
		y_first_nn = y_temp(idx,col)
		d_fnn(a,col) = sqrt(((x_temp(a,col)-x_first_nn)**2)+ &
			((y_temp(a,col)- y_first_nn)**2));
		!write(*,*) d_fnn(a,col)
		a = a + 1
	end do
		!write(*,*) d_fnn(1:N,col) 

	!calculating delta_phi----------------------------------------
	a = 1
	do while (a .le. N)
		if (a .ne. idx_new) then
			phi_other(a,col) = atan2(y_temp(a,col)-	     &
			y_temp(a,col-1),x_temp(a,col)-x_temp(a,col-1))
		else
			phi_focal(a,col) = atan2(y_temp(a,col)-	     &
			y_temp(a,col-1),x_temp(a,col)-x_temp(a,col-1))
		end if 
		a = a + 1
	end do
	!write(*,*) phi_focal(:,col),phi_other(:,col)

	a = 1
	do while (a .le. N)
		if (a .ne. idx_new) then
			delta_phi(a,col) = phi_other(a,col)-	     &
					phi_focal(idx_new,col)
			do while (delta_phi(a,col) .gt. pi)
				delta_phi(a,col) = delta_phi(a,col)- &
							(2.0*pi)
			end do
			do while (delta_phi(a,col) .lt. -pi)
				delta_phi(a,col) = delta_phi(a,col)+ &
							(2.0*pi)
			end do
		delta_phi(a,col) = -1.0*delta_phi(a,col)
		end if
		!write (*,*) delta_phi(a,col)
		
		a = a + 1 
	end do
	!write (*,*) delta_phi(:,col)
	!*************************************************************

	!calculating psi----------------------------------------------
	a = 1
	do while (a .le. N)
		if (a .ne. idx_new) then
			psi(a,col) = atan2(y_temp(a,col)-	     &
					y_temp(idx_new,col),	     &
				(x_temp(a,col)-x_temp(idx_new,col))) &
				- phi_focal(idx_new,col)
		do while (psi(a,col) .gt. pi) 
			psi(a,col) = psi(a,col) - (2.0*pi)
		end do
		do while (psi(a,col) .lt. -pi)
			psi(a,col) = psi(a,col) + (2.0*pi)
		end do
		psi(a,col) = -1.0*psi(a,col)
		end if
		a = a + 1
		!write(*,*) psi(a,col)
	end do
	!write(*,*) psi(:,col)
	!*************************************************************

	!substituting the angles obtained to get the updated angle
	g(idx_new,col) = gaus_rand()
	!write(*,*) gaus_rand()
	a = 1
	del_phi_ali_tot = 0.0
	del_phi_att_tot = 0.0
	do while (a .le. N)
		!write (*,*) gamma_att
		!write (*,*) gamma_ali
		if (a .ne. idx_new) then
			!write(*,*) psi(a,col)
			del_phi_att_tot = del_phi_att_tot +	     &
	del_phi_att(gamma_att(e),d(a,col),psi(a,col),delta_phi(a,col))
			do while (del_phi_att_tot .gt. pi)
				del_phi_att_tot = del_phi_att_tot    &
					- (2.0*pi)
			end do
			do while (del_phi_att_tot .lt. -pi)
				del_phi_att_tot = del_phi_att_tot    &
					+ (2.0*pi)
			end do
		end if
		if (a .ne. idx_new) then
			del_phi_ali_tot = del_phi_ali_tot +	     &
	del_phi_ali(gamma_ali,d(a,col),psi(a,col),delta_phi(a,col))
			do while (del_phi_ali_tot .gt. pi)
				del_phi_ali_tot = del_phi_ali_tot    &
					- (2.0*pi)
			end do
			do while (del_phi_ali_tot .lt. -pi)
				del_phi_ali_tot = del_phi_ali_tot    &
					+ (2.0*pi)
			end do
		end if 
		a = a + 1
	end do
	!write(*,*) del_phi_att_tot
	!write(*,*) del_phi_ali_tot
	phi(idx_new,col) = phi(idx_new,col-1) + del_phi_att_tot +    &
			del_phi_ali_tot + (gamma_rand(f)*g(idx_new,col))
	do while (phi(idx_new,col) .gt. pi)
		phi(idx_new,col) = phi(idx_new,col) - (2.0*pi)
	end do
	do while (phi(idx_new,col) .lt. -pi)
		phi(idx_new,col) = phi(idx_new,col) + (2.0*pi)
	end do
	!*************************************************************

	!preparing the arrays for the next loop----------------------- 
	a = 1 
	do while (a .le. N)
		if (a .ne. idx_new) then
			v(a,col) = v(a,col-1)
			tau(a,col) = tau(a,col-1)
			phi(a,col) = phi(a,col-1)
			time(a,col+1) = time(a,col)
		else
			call random_number(r)
			v(a,col) = generate_v(r)
			call random_number(r)
			!write(*,*) r
			tau(a,col) = generate_tau(r)
			time(a,col+1) = time(a,col) + tau(a,col)
		end if
		a = a + 1
	end do
	!*************************************************************
	
	!to find polarisation-----------------------------------------
	pol_sum = 0.0
	cos_term = 0.0
	sin_term = 0.0
	do a = 1,N,1
		cos_term = cos_term + cos(phi(a,col))
		sin_term = sin_term + sin(phi(a,col))
	end do
	!write(*,*) cos_term,sin_term
	pol_mod = sqrt(((cos_term)**2)+((sin_term)**2))
	pol(col) = (1.0/N)*pol_mod
	!write(*,*) pol(col),(1.0/N)*pol_mod,pol_mod
	!*************************************************************

	!to find the radius of gyration-------------------------------
	summ = 0.0
	do a=1,N,1
		do c=a+1,N,1
			summ= summ+(((x_temp(a,col)-x_temp(c,col))**2)&
		+((y_temp(a,col)-y_temp(c,col))**2))
		end do
	end do
	r_sq(col) = (2.0/(N*(N-1)))*summ
	write(13,*) time_mod(1,col),r_sq(col)
	!*************************************************************
	
	col = col + 1
	!write(*,*) col
end do
	columns(i) = col-1
	write(*,*) col-1
	close(13)
end do !end for each sample
	open(unit=9,file="col_minval.txt",status="replace",	     &
		action="write",position="rewind",		     &
iostat=column_file_status)
	write(*,*) "column_file_status=", column_file_status
	write(9,*) minval(columns)
	close(9)

call chdir('..',status = chdir_status)
write(*,*) "chdir_status=",chdir_status
end do !end for gamma_att
end do !end for gamma_rand

!storing the trajectories---------------------------------------------
!do i = 1,col-1
!	do j=1,N
!		if (j .lt. N) then
!			write(10,'(f10.5,7x)',advance="no")    &
!				x_temp(j,i)
!			write(10,'(f10.5,7x)',advance="no")    &
!				y_temp(j,i)
!		else
!			write(10,'(f10.5,7x)',advance="no")    &
!				x_temp(j,i)
!			write(10,'(f10.5)')		       &
!				y_temp(j,i)
!		end if
!	end do
!end do
!*********************************************************************

!storing dist_fnn-----------------------------------------------------
!summ = 0
!do i=1,col-1
!	summ = summ + (sum(d_fnn(1:N,i))/N)
!	!write(*,*) sum(d_fnn(1:N,i))/N
!end do
!write(*,*) "fnn=",summ/(col-1)
!*********************************************************************

!storing polarsation--------------------------------------------------
!summ=0
!do i=1,col-1
!	summ = summ + pol(i)
!end do
!write(*,*) summ/(col-1)
!*********************************************************************

!storing sqrt(<r^2>/N)-----------------------------------------------
!summ = 0
!do i = 1,col-1
!	summ = summ + r_sq(i)
!end do
!write(*,*) sqrt(summ/(N*(col-1)))
!*********************************************************************

!close(10)
!close(11)
!close(12)
end program nfish
!*********************************************************************

!*********************************************************************
!RANDOM SEED GENERATOR FOR THE RANDOM_NUMBER PROGRAM
!*********************************************************************

subroutine init_random_seed()
	use iso_fortran_env, only: int64
	implicit none
	integer, allocatable :: seed(:)
	integer :: i, n, un, istat, dt(8), pid
	integer(int64) :: t
		
	call random_seed(size = n)
	allocate(seed(n))
	! First try if the OS provides a random number generator
	open(newunit=un, file="/dev/urandom", access="stream",        &
	form="unformatted", action="read", status="old", iostat=istat)
	if (istat == 0) then
		read(un) seed
		close(un)
	else
		! Fallback to XOR:ing the current time and pid. The PID is
		! useful in case one launches multiple instances of the same
		! program in parallel.
		call system_clock(t)
		if (t == 0) then
			call date_and_time(values=dt)
			t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 &
							* 1000        &
			+ dt(2) * 31_int64 * 24 * 60 * 60 * 1000      &
			+ dt(3) * 24_int64 * 60 * 60 * 1000           &
			+ dt(5) * 60 * 60 * 1000                      &
			+ dt(6) * 60 * 1000 + dt(7) * 1000            &
			+ dt(8)
		end if
		pid = getpid()
		t = ieor(t, int(pid, kind(t)))
		do i = 1, n
			seed(i) = lcg(t)
		end do
		end if
		call random_seed(put=seed)	!value of 'put' is the array used to initialize this thread's seed.
	contains
		! This simple PRNG might not be good enough for real work, but is
		! sufficient for seeding a better PRNG.
		function lcg(s)
			integer :: lcg
			integer(int64) :: s
			if (s == 0) then
				s = 104729
			else
				s = mod(s, 4294967296_int64)
			end if
				s = mod(s * 279470273_int64,          &
				4294967291_int64)
				lcg = int(mod(s, int(huge(0),         &
 				int64)), kind(0))
		end function lcg
end subroutine init_random_seed
!*********************************************************************
