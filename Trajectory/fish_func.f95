!*********************************************************************
!Written by Rachith Aiyappa at IISc, Bangalore
!Email rachithaiyappa96@gmail.com for any queries
!This module contains the functions for the N fish simulation 
!*********************************************************************
!
module fish_func
contains
!*********************************************************************
!function generate_v
!This function generates kicking speeds of the respective focal fish
!Units of BL_tau
!*********************************************************************
function generate_v(r)
	real,parameter :: pi = 3.14159265 
	real :: r
	generate_v = sqrt((-16/pi)*log(1-r))
end function generate_v
!---------------------------------------------------------------------
!
!********************************************************************!
!function generate_tau
!This function generates kicking times of the respective focal fish
!Units of BL_tau
!********************************************************************!
function generate_tau(r)
	real,parameter :: pi = 3.14159265
	real :: r	
	generate_tau = sqrt((-4/pi)*log(1-r))
end function generate_tau
!--------------------------------------------------------------------
!
!********************************************************************
!function dist_calc
!This function calculates the distance from the focal fish
!********************************************************************
function dist_calc(col,focal_x,focal_y,x_temp,y_temp,N) result(d)
	integer,intent(in) :: col,N
	integer :: i
	
	real, dimension(N,1) :: d
	real,intent(in),dimension(N,1) :: x_temp,y_temp
	real,intent(in) :: focal_x,focal_y
	
	!write (*,*) focal_x,x_temp
	do i=1,N,1
		!write(*,*) d(i,col)
		d(i,1) = sqrt(((focal_x-x_temp(i,1))**2)	&
				+ ((focal_y-y_temp(i,1))**2))
		!write(*,*) d(i,col)
	end do
	!write(*,*) d
end function
!--------------------------------------------------------------------
!
!********************************************************************
!function del_phi_att
!This function calculate the change in heading angle due to attraction
!********************************************************************
function del_phi_att(gamma_att,d,psi,delta_phi)
	!a is included to restore anticlockwise angle as positive
	real,parameter :: a=-1.0
	real,parameter :: l_att=5.0,d_zero=1.0
	!real,parameter :: l_att=1.0 !proposed by Clement
	real :: f_att,o_att,e_att
	real,intent(in) :: gamma_att,d,psi,delta_phi
	
	!write(*,*) gamma_att, d, psi, delta_phi
	!write(*,*) 

	!original functions in the paper
	f_att = ((d-d_zero)/d_zero)/(1.0+((d/l_att)**2.0))
	o_att = sin(a*psi)*(1-(0.33*cos(a*psi)))
	e_att = 1.0-(0.48*cos(a*delta_phi))-			    &
				(0.31*cos(2.0*(a*delta_phi)))
	

	!functions proposed by Clement
	!f_att = (d/d_zero)/(1.0+((d/l_att)**2.0))
	!o_att = sin(a*psi)*(1+cos(a*psi))
	!e_att = 1
	
	del_phi_att = gamma_att*f_att*o_att*e_att
end function del_phi_att
!---------------------------------------------------------------------

!********************************************************************
!function del_phi_ali
!This function calculate the change in heading angle due to alignment
!********************************************************************
function del_phi_ali(gamma_ali,d,psi,delta_phi)
	
	!a is included to restore anticlockwise angle as positive
	real,parameter :: a=-1.0
	real,parameter :: l_ali = 5.0
	real :: f_ali,o_ali,e_ali
	real,intent(in) :: gamma_ali,d,psi,delta_phi
	
	!write(*,*) gamma_ali, d, psi, delta_phi

	!original funcions in the paper
	f_ali = exp(-1.0*((d/l_ali)**2.0))
	e_ali = (1.0) + (0.6*cos(a*psi))-(0.32*cos(2.0*(a*psi)))
	o_ali = (sin(a*delta_phi))*((1.0)+			    &
			(0.30*cos(2.0*(a*delta_phi))))

	!functions proposed by Clement
	!f_ali = exp(-1.0*((d/l_ali)**2.0))
	!o_ali = sin(a*delta_phi)
	!e_ali = 1+cos(a*psi)

	del_phi_ali = gamma_ali*f_ali*o_ali*e_ali
end function del_phi_ali
!---------------------------------------------------------------------

!********************************************************************
!function gaus_rand
!This uses the marsaglia polar method to create (pseudo-optional) 
!random number pairs that have a normal distribution. The function
!saves one of the random numbers from the generated pair as a spare 
!to be returned from the generated pair for the next call to the func
!Returns a real scalar
!********************************************************************
function gaus_rand()
	real :: gaus_rand
	real :: mean=0,sd=1
	!real,intent(inout) :: x,y
	real :: r,x,y
	real,save :: spare
	logical, save :: has_spare
	call init_random_seed()
	!use a spare saved from a previous run if it exits
	if(has_spare) then
		has_spare = .false.
		gaus_rand = mean + (sd*spare)
		return
	else
		r = 1.0
		do while (r .ge. 1.0)
			!spreading distribution to [-1,1]
			call random_number(x)
			call random_number(y)
			x = (x*2.0)-1.0
			y = (y*2.0)-1.0
			r = (x*x) + (y*y)
		end do

		r = sqrt((-2.0*log(r))/r)
		gaus_rand = mean + (sd*x*r)
		spare = y * r
		has_spare = .true.
		return
	end if
end function
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

end module fish_func
