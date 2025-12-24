program SW1L
	implicit none
	
	real(8) :: pi, omega, g, Lx, Ly, dx, dy, dt, time, f
	integer :: i, j, k, Nx, Ny, Nt, tspawn

	character(50) :: IC, scheme
	logical :: rotating_frame 

	real(8), allocatable :: x(:), y(:), t(:), u(:,:,:), v(:,:,:), eta(:,:,:), H(:,:), pertur(:,:)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	pi = 3.14159265359
	omega = 2*pi/(24*3600)	
	g = 9.81
	tspawn = 2

	IC = 'detroit' ! choose : 'center', 'island', 'detroit'
	rotating_frame = .true.
	scheme = 'LPEF' ! chosse : 'EF', 'LP', 'LPEF'


	Lx = 5000
	Ly = 5000

	Nx = 50
	Ny = 50

	dx = Lx/Nx
	dy = Ly/Ny

	dt = 20
	time = 5000
	!time = 7000
	Nt = time/dt

	allocate(t(Nt), x(Nx), y(Ny))

	t(:) = 0
	x(:) = 0
	y(:) = 0

	do i=1,Nt-1
		t(i+1) = t(i) + dt
	end do
	do j=1,Nx-1
		x(j+1) = x(j) + dx
	end do
	do k=1,Ny-1
                y(k+1) = y(k) + dy
        end do

	allocate(u(Nt, Nx, Ny), v(Nt, Nx, Ny), eta(Nt, Nx, Ny), H(Nx, Ny), pertur(Nx, Ny))

	H(:,:) = 1

	! initial perturbation

	call gaussian2D(x, y, eta, pi, Nx, Ny, pertur, IC)
	eta(1,:,:) = pertur*5*1e5
	
	!u(1,:10,:) = 0.25

	
	if (rotating_frame .eqv. .false.) then
		f = 0.
		print *, 'no coriolis'
		
		if (scheme == 'EF') then
			call integrator_EF(Nt, Nx, Ny, u, v, eta, dx, dy, dt, g, H, f, IC)
		elseif (scheme == 'LP') then
			call integrator_LP(Nt, Nx, Ny, u, v, eta, dx, dy, dt, g, H, f, IC)
		elseif (scheme == 'LPEF') then		
			call integrator_LP_EF(Nt, Nx, Ny, u, v, eta, dx, dy, dt, g, H, f, IC)
		endif


	elseif (rotating_frame .eqv. .true.) then
		print *, Nt
		f = 2*omega*sin(45.)
		print *, f
		
		if (scheme == 'EF') then
                        call integrator_EF(Nt, Nx, Ny, u, v, eta, dx, dy, dt, g, H, f, IC)
                elseif (scheme == 'LP') then
                        call integrator_LP(Nt, Nx, Ny, u, v, eta, dx, dy, dt, g, H, f, IC)
		elseif (scheme == 'LPEF') then
			call integrator_LP_EF(Nt, Nx, Ny, u, v, eta, dx, dy, dt, g, H, f, IC)

                endif
	end if


	print *, eta


contains

	subroutine gaussian2D(x, y, eta, pi, Nx, Ny, pertur, IC)

		integer :: j, k
		real(8), allocatable :: xx(:,:), yy(:,:)		
		real(8) :: sx, sy, xbar, ybar, Vx, Vy
		
		integer, intent(in) :: Nx, Ny
		character(50), intent(in) :: IC
		real(8), intent(inout) :: x(:), y(:), eta(:,:,:), pi

		real(8), allocatable, intent(out) :: pertur(:,:)


		! meshgrid
		allocate(pertur(Nx,Ny))

		! variances

		if (IC == 'center') then
			xbar = sum(x)/(Nx)
			ybar = sum(y)/(Ny)
		elseif (IC == 'island') then
			xbar = sum(x)/(Nx/1.8)    
                        ybar = sum(y)/(Ny/1.8)
		elseif (IC == 'detroit') then
			xbar = sum(x)/(Nx/2.2)
			ybar = sum(y)/Ny
		end if



		Vx = sum((x - xbar)**2)/Nx
		Vy = sum((y - ybar)**2)/Ny

		!std
		sx = 0.25*(Vx**(0.5)) 
		sy = 0.25*(Vy**(0.5)) 

		! pas le droit en f90 de rentrer une valeur puis de la modifier comme ça
		! stocker dans un tampon
		
		do j=1,Nx
			do k=1,Ny
				!pertur(j,k) = eta(1,j,k) + (1/(2*pi*sx*sy))*exp(-0.5*(( (x(j)-xbar)/sx )**2 + ((y(k)-ybar)/sy)**2 ))
				pertur(j,k) = eta(1,j,k) + (1/(2*pi*sx*sy))*exp(-0.5*(( (x(j)-xbar)/sx )**2 + ((y(k)-ybar)/sy)**2 ))
			end do
		end do


	end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EULER

	subroutine integrator_EF(Nt, Nx, Ny, u, v, eta, dx, dy, dt, g, H, f, IC)

		integer :: i, j, k

                integer, intent(in) :: Nt, Nx, Ny
                real(8), intent(in) :: dx, dy, dt, g, f
                character(50), intent(in) :: IC

                real(8), intent(inout) :: u(:,:,:), v(:,:,:), eta(:,:,:), H(:,:)

                do i=1,Nt-1
                        do j=2,Nx-1
                                do k=2,Ny-1

                                        !print *,'Euler'
					u(i+1,j,k) = u(i,j,k)-dt*(g* (eta(i,j+1,k)-eta(i,j-1,k))/(2*dx) + f*v(i,j,k))
					v(i+1,j,k) = v(i,j,k)-dt*(g* (eta(i,j,k+1)-eta(i,j,k-1))/(2*dy) - f*u(i,j,k))
					eta(i+1,j,k) = eta(i,j,k)-dt*(H(j,k)*((u(i,j+1,k)-u(i,j-1,k))/(2*dx) + (v(i,j,k+1)-v(i,j,k-1))/(2*dy)))

				end do
                        end do

                        ! BC's : rigid wall

                        u(i,1,:) = 0.
                        u(i,Nx,:) = 0.
                        v(i,:,1) = 0.
                        v(i,:,Ny) = 0.
			
			u(i+1,1,:)  = 0.
                        u(i+1,Nx,:) = 0.
                        v(i+1,:,1)  = 0.
                        v(i+1,:,Ny) = 0.

                        ! simple neuman
                        eta(i+1,1,:)  = eta(i+1,2,:)
                        eta(i+1,Nx,:) = eta(i+1,Nx-1,:)
                        eta(i+1,:,1)  = eta(i+1,:,2)
                        eta(i+1,:,Ny) = eta(i+1,:,Ny-1)


                        if (IC == 'island') then
                                ! island
                                u(i,Ny/3+5:Ny/3+15,20:30) = 0.
                                v(i,Ny/3+5:Ny/3+15,20:30) = 0.

                                u(i+1,Ny/3+5:Ny/3+15,20:30) = 0
                                v(i+1,Ny/3+5:Ny/3+15,20:30) = 0
                        elseif (IC == 'detroit') then
                                ! COTE GAUCHE
                                u(i,Ny/3+5:Ny/3+15,1:24) = 0
                                v(i,Ny/3+5:Ny/3+15,1:24) = 0

                                u(i+1,Ny/3+5:Ny/3+15,1:24) = 0
                                v(i+1,Ny/3+5:Ny/3+15,1:24) = 0

                                ! COTE DROITE
                                u(i,Ny/3+5:Ny/3+15,26:50) = 0
                                v(i,Ny/3+5:Ny/3+15,26:50) = 0

                                u(i+1,Ny/3+5:Ny/3+15,26:50) = 0
                                v(i+1,Ny/3+5:Ny/3+15,26:50) = 0

                        end if


                end do

	end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LEAP-FROG



	subroutine integrator_LP(Nt, Nx, Ny, u, v, eta, dx, dy, dt, g, H, f, IC)

		integer :: tspawn, i, j, k

		integer, intent(in) :: Nt, Nx, Ny
		real(8), intent(in) :: dx, dy, dt, g, f
		character(50), intent(in) :: IC
		
		real(8), intent(inout) :: u(:,:,:), v(:,:,:), eta(:,:,:), H(:,:)

		tspawn = 3

		do i=1,Nt-1					
			if (i <= tspawn) then
				do j=2,Nx-1
				!print *, 'Euler init'
					do k=2,Ny-1
						u(i+1,j,k) = u(i,j,k)-dt*(g* (eta(i,j+1,k)-eta(i,j-1,k))/(2*dx) + f*v(i,j,k))
						v(i+1,j,k) = v(i,j,k)-dt*(g* (eta(i,j,k+1)-eta(i,j,k-1))/(2*dy) - f*u(i,j,k))
						eta(i+1,j,k) = eta(i,j,k)-dt*(H(j,k)*((u(i,j+1,k)-u(i,j-1,k))/(2*dx) + (v(i,j,k+1)-v(i,j,k-1))/(2*dy)))
					end do
				end do

			else 
				do j=2,Nx-1
					!print *, 'Leap-Frog'
					do k=2,Ny-1
						u(i+1,j,k) = u(i-1,j,k)-2*dt*(g* (eta(i,j+1,k)-eta(i,j-1,k))/(2*dx) + f*v(i,j,k))
						v(i+1,j,k) = v(i-1,j,k)-2*dt*(g* (eta(i,j,k+1)-eta(i,j,k-1))/(2*dy) - f*u(i,j,k))
						eta(i+1,j,k) = eta(i-1,j,k)-2*dt*(H(j,k)*((u(i,j+1,k)-u(i,j-1,k))/(2*dx) + (v(i,j,k+1)-v(i,j,k-1))/(2*dy)))
					end do
				end do

			end if


			! BC's : rigid wall

			u(i,1,:) = 0.
			u(i,Nx,:) = 0.
			v(i,:,1) = 0.
			v(i,:,Ny) = 0.

			u(i+1,1,:)  = 0.
			u(i+1,Nx,:) = 0.
			v(i+1,:,1)  = 0.
			v(i+1,:,Ny) = 0.

			! simple neuman
			eta(i+1,1,:)  = eta(i+1,2,:)
			eta(i+1,Nx,:) = eta(i+1,Nx-1,:)
			eta(i+1,:,1)  = eta(i+1,:,2)
			eta(i+1,:,Ny) = eta(i+1,:,Ny-1)


			if (IC == 'island') then
				! island
				u(i,Ny/3+5:Ny/3+15,20:30) = 0.
				v(i,Ny/3+5:Ny/3+15,20:30) = 0.

				u(i+1,Ny/3+5:Ny/3+15,20:30) = 0
				v(i+1,Ny/3+5:Ny/3+15,20:30) = 0
			elseif (IC == 'detroit') then
				! COTE GAUCHE
				u(i,Ny/3+5:Ny/3+15,1:24) = 0
				v(i,Ny/3+5:Ny/3+15,1:24) = 0

				u(i+1,Ny/3+5:Ny/3+15,1:24) = 0
                                v(i+1,Ny/3+5:Ny/3+15,1:24) = 0

				! COTE DROITE
				u(i,Ny/3+5:Ny/3+15,26:50) = 0
				v(i,Ny/3+5:Ny/3+15,26:50) = 0

				u(i+1,Ny/3+5:Ny/3+15,26:50) = 0
                                v(i+1,Ny/3+5:Ny/3+15,26:50) = 0

			end if


		end do
	
	end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MIXTE : LEAP-FROG EULER




subroutine integrator_LP_EF(Nt, Nx, Ny, u, v, eta, dx, dy, dt, g, H, f, IC)

                integer :: tspawn, i, j, k
		integer :: ct, inter, r    ! integers for the mix scheme


                integer, intent(in) :: Nt, Nx, Ny
                real(8), intent(in) :: dx, dy, dt, g, f
                character(50), intent(in) :: IC

                real(8), intent(inout) :: u(:,:,:), v(:,:,:), eta(:,:,:), H(:,:)

                tspawn = 3
		inter = 40
		ct = 1
		r = 0



                do i=1,Nt-1
                        
                        if (i <= tspawn) then
				do j=2,Nx-1
					!print *, 'Euler init'
					do k=2,Ny-1
						u(i+1,j,k) = u(i,j,k)-dt*(g* &
							(eta(i,j+1,k)-eta(i,j-1,k))/(2*dx) + f*v(i,j,k))
						v(i+1,j,k) = v(i,j,k)-dt*(g* &
							(eta(i,j,k+1)-eta(i,j,k-1))/(2*dy) - f*u(i,j,k))
						eta(i+1,j,k) = eta(i,j,k)-dt*(H(j,k)* &
							((u(i,j+1,k)-u(i,j-1,k))/(2*dx) + (v(i,j,k+1)-v(i,j,k-1))/(2*dy)))

					end do
				end do

                        else
				if (i==inter*ct) then
					do while (r < 10)
						do j=2,Nx-1
							!print *, 'Euler'
							do k=2,Ny-1

								u(i+1,j,k) = u(i,j,k)-dt*(g* &
									(eta(i,j+1,k)-eta(i,j-1,k))/(2*dx) + f*v(i,j,k))
								v(i+1,j,k) = v(i,j,k)-dt*(g* &
									(eta(i,j,k+1)-eta(i,j,k-1))/(2*dy) - f*u(i,j,k))
								eta(i+1,j,k) = eta(i,j,k)-dt*(H(j,k)* &
									((u(i,j+1,k)-u(i,j-1,k))/(2*dx) + (v(i,j,k+1)-v(i,j,k-1))/(2*dy)))
							end do
						end do

						r = r+1


					end do

					ct = ct + 1
					r = 0

				else

					do j=2,Nx-1
						!print *, 'Leap-Frog'
						do k=2,Ny-1

							u(i+1,j,k) = u(i-1,j,k)-2*dt*(g* &
								(eta(i,j+1,k)-eta(i,j-1,k))/(2*dx) + f*v(i,j,k))
							v(i+1,j,k) = v(i-1,j,k)-2*dt*(g* &
								(eta(i,j,k+1)-eta(i,j,k-1))/(2*dy) - f*u(i,j,k))
							eta(i+1,j,k) = eta(i-1,j,k)-2*dt*(H(j,k)* &
								((u(i,j+1,k)-u(i,j-1,k))/(2*dx) + (v(i,j,k+1)-v(i,j,k-1))/(2*dy)))
                        
						end do
					end do
				end if

			end if

        
                        ! BC's : rigid wall

                        u(i,1,:) = 0.
                        u(i,Nx,:) = 0.
                        v(i,:,1) = 0.
                        v(i,:,Ny) = 0.

			u(i+1,1,:)  = 0.
                        u(i+1,Nx,:) = 0.
                        v(i+1,:,1)  = 0.
                        v(i+1,:,Ny) = 0.

                        ! simple neuman
                        eta(i+1,1,:)  = eta(i+1,2,:)
                        eta(i+1,Nx,:) = eta(i+1,Nx-1,:)
                        eta(i+1,:,1)  = eta(i+1,:,2)
                        eta(i+1,:,Ny) = eta(i+1,:,Ny-1)


                        if (IC == 'island') then
                                ! island
                                u(i,Ny/3+5:Ny/3+15,20:30) = 0.
                                v(i,Ny/3+5:Ny/3+15,20:30) = 0.

                                u(i+1,Ny/3+5:Ny/3+15,20:30) = 0
                                v(i+1,Ny/3+5:Ny/3+15,20:30) = 0
                        elseif (IC == 'detroit') then
                                ! COTE GAUCHE
                                u(i,Ny/3+5:Ny/3+15,1:24) = 0
                                v(i,Ny/3+5:Ny/3+15,1:24) = 0

                                u(i+1,Ny/3+5:Ny/3+15,1:24) = 0
                                v(i+1,Ny/3+5:Ny/3+15,1:24) = 0

                                ! COTE DROITE
                                u(i,Ny/3+5:Ny/3+15,26:50) = 0
                                v(i,Ny/3+5:Ny/3+15,26:50) = 0

                                u(i+1,Ny/3+5:Ny/3+15,26:50) = 0
                                v(i+1,Ny/3+5:Ny/3+15,26:50) = 0

                        end if


                end do

	end subroutine








end program
