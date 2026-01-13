program SW2L
	implicit none
	
	real(8) :: pi, omega, g, g2, Lx, Ly, dx, dy, dt, time, f, rho1, rho2, rho0
	integer :: i, j, k, Nx, Ny, Nt, tspawn

	character(50) :: IC, scheme
	logical :: rotating_frame 

	real(8), allocatable :: x(:), y(:), t(:), u1(:,:,:), v1(:,:,:), eta1(:,:,:), H1(:,:), pertur(:,:)
	real(8), allocatable :: u2(:,:,:), v2(:,:,:), eta2(:,:,:), H2(:,:)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	pi = 3.14159265359
	omega = 2*pi/(24*3600)	
	g = 9.81
	tspawn = 2

	rho2 = 1020
	rho1 = 1010
	rho0 = 1013
	g2= g*(rho1 - rho2)/rho0  !reduced gravity


	IC = 'center' ! choose : 'center', 'island', 'detroit'
	rotating_frame = .true.
	scheme = 'LP' ! chosse : 'EF', 'LP', 'LPEF'


	Lx = 5000
	Ly = 5000

	Nx = 50
	Ny = 50

	dx = Lx/(Nx-1)
	dy = Ly/(Ny-1)

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

	allocate(u1(Nt, Nx, Ny), v1(Nt, Nx, Ny), eta1(Nt, Nx, Ny), H1(Nx, Ny), pertur(Nx, Ny))
	allocate(u2(Nt, Nx, Ny), v2(Nt, Nx, Ny), eta2(Nt, Nx, Ny), H2(Nx, Ny))

	H1(:,:) = 1.
	H2(:,:) = 1.

	! initial perturbation

	call gaussian2D(x, y, eta1, pi, Nx, Ny, pertur, IC)
	!eta1(1,:,:) = pertur*5*1e5
	eta2(1,:,:) = -pertur*1e5	


	!u(1,:10,:) = 0.25

	
	if (rotating_frame .eqv. .false.) then
		f = 0.
		print *, 'no coriolis'
		call integrator_LP(Nt, Nx, Ny, u1, u2, v1, v2, eta1, eta2, dx, dy, dt, g, H1, H2, f, IC)
		
	elseif (rotating_frame .eqv. .true.) then
		print *, Nt
		f = 2*omega*sin(45.)
		print *, f
		
		call integrator_LP(Nt, Nx, Ny, u1, u2, v1, v2, eta1, eta2, dx, dy, dt, g, H1, H2, f, IC)
	end if


	print *, eta1
	print *, eta2


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
! LEAP-FROG



	subroutine integrator_LP(Nt, Nx, Ny, u1, u2, v1, v2, eta1, eta2, dx, dy, dt, g, H1, H2, f, IC)

		integer :: tspawn, i, j, k

		integer, intent(in) :: Nt, Nx, Ny
		real(8), intent(in) :: dx, dy, dt, g, f
		character(50), intent(in) :: IC
		
		real(8), intent(inout) :: u1(:,:,:), v1(:,:,:), eta1(:,:,:), H1(:,:)
		real(8), intent(inout) :: u2(:,:,:), v2(:,:,:), eta2(:,:,:), H2(:,:)
	
		real(8), allocatable :: detadt(:,:,:)

		allocate(detadt(Nt,Nx,Ny))


		tspawn = 3

		do i=1,Nt-1					
			if (i <= tspawn) then
				do j=2,Nx-1
				!print *, 'Euler init'
					do k=2,Ny-1
						u1(i+1,j,k)=u1(i,j,k)-dt*(g* (eta1(i,j+1,k)-eta1(i,j-1,k))/(2*dx) + f*v1(i,j,k))
						v1(i+1,j,k)=v1(i,j,k)-dt*(g* (eta1(i,j,k+1)-eta1(i,j,k-1))/(2*dy) - f*u1(i,j,k))
						detadt(i,j,k)=-H2(j,k)*((u2(i,j+1,k) - u2(i,j-1,k))/(2*dx) + (v2(i,j,k+1) - v2(i,j,k-1))/(2*dy))
						eta1(i+1,j,k)=eta1(i,j,k)-dt*(H1(j,k)*(u1(i,j+1,k)-u1(i,j-1,k))/(2*dx) + (v1(i,j,k+1)-v1(i,j,k-1))/(2*dy)) - detadt(i,j,k)
                        
						! 2nd layer
						u2(i+1,j,k)=u2(i,j,k)-dt*(g* (eta1(i,j+1,k)-eta1(i,j-1,k))/(2*dx) - g2*(eta2(i,j+1,k)-eta2(i,j-1,k))/(2*dx)+f*v2(i,j,k))
						v2(i+1,j,k)=v2(i,j,k)-dt*(g* (eta1(i,j,k+1)-eta1(i,j,k+1))/(2*dy) - g2*(eta2(i,j,k+1)-eta2(i,j,k+1))/(2*dy)-f*u2(i,j,k))
						eta2(i+1,j,k)=eta2(i,j,k)-dt*H2(j,k)*((u2(i,j+1,k)-u2(i,j-1,k))/(2*dx) + (v2(i,j,k+1)-v2(i,j,k-1))/(2*dy))

					end do
				end do

			else 
				do j=2,Nx-1
					!print *, 'Leap-Frog'
					do k=2,Ny-1
						u1(i+1,j,k)=u1(i-1,j,k)-2*dt*(g* (eta1(i,j+1,k)-eta1(i,j-1,k))/(2*dx) + f*v1(i,j,k))
						v1(i+1,j,k)=v1(i-1,j,k)-2*dt*(g* (eta1(i,j,k+1)-eta1(i,j,k-1))/(2*dy) - f*u1(i,j,k))
						detadt(i,j,k)=-H2(j,k)*((u2(i,j+1,k) - u2(i,j-1,k))/(2*dx) + (v2(i,j,k+1) - v2(i,j,k-1))/(2*dy))
						eta1(i+1,j,k)=eta1(i-1,j,k)-2*dt*(H1(j,k)*(u1(i,j+1,k)-u1(i,j-1,k))/(2*dx) + (v1(i,j,k+1)-v1(i,j,k-1))/(2*dy)) - detadt(i,j,k)
                        
						! 2nd layer
						u2(i+1,j,k)=u2(i-1,j,k)-2*dt*(g* (eta1(i,j+1,k)-eta1(i,j-1,k))/(2*dx) - g2*(eta2(i,j+1,k)-eta2(i,j-1,k))/(2*dx)+f*v2(i,j,k))
						v2(i+1,j,k)=v2(i-1,j,k)-2*dt*(g* (eta1(i,j,k+1)-eta1(i,j,k+1))/(2*dy) - g2*(eta2(i,j,k+1)-eta2(i,j,k+1))/(2*dy)-f*u2(i,j,k))
						eta2(i+1,j,k)=eta2(i-1,j,k)-2*dt*H2(j,k)*((u2(i,j+1,k)-u2(i,j-1,k))/(2*dx) + (v2(i,j,k+1)-v2(i,j,k-1))/(2*dy))
                       
					end do
				end do

			end if


			! BC's : rigid wall

			u1(i,1,:) = 0.
			u1(i,Nx,:) = 0.
			v1(i,:,1) = 0.
			v1(i,:,Ny) = 0.

			u1(i+1,1,:)  = 0.
			u1(i+1,Nx,:) = 0.
			v1(i+1,:,1)  = 0.
			v1(i+1,:,Ny) = 0.

			u2(i,1,:) = 0.
                        u2(i,Nx,:) = 0.
                        v2(i,:,1) = 0.
                        v2(i,:,Ny) = 0.

                        u2(i+1,1,:)  = 0.
                        u2(i+1,Nx,:) = 0.
                        v2(i+1,:,1)  = 0.
                        v2(i+1,:,Ny) = 0.

			! simple neuman
			eta1(i+1,1,:)  = eta1(i+1,2,:)
			eta1(i+1,Nx,:) = eta1(i+1,Nx-1,:)
			eta1(i+1,:,1)  = eta1(i+1,:,2)
			eta1(i+1,:,Ny) = eta1(i+1,:,Ny-1)

			eta2(i+1,1,:)  = eta2(i+1,2,:)
                        eta2(i+1,Nx,:) = eta2(i+1,Nx-1,:)
                        eta2(i+1,:,1)  = eta2(i+1,:,2)
                        eta2(i+1,:,Ny) = eta2(i+1,:,Ny-1)


			if (IC == 'island') then
				! island
				u1(i,Ny/3+5:Ny/3+15,20:30) = 0.
				v1(i,Ny/3+5:Ny/3+15,20:30) = 0.

				u1(i+1,Ny/3+5:Ny/3+15,20:30) = 0
				v1(i+1,Ny/3+5:Ny/3+15,20:30) = 0

				u2(i,Ny/3+5:Ny/3+15,20:30) = 0.
                                v2(i,Ny/3+5:Ny/3+15,20:30) = 0.

                                u2(i+1,Ny/3+5:Ny/3+15,20:30) = 0
                                v2(i+1,Ny/3+5:Ny/3+15,20:30) = 0


			elseif (IC == 'detroit') then
				! COTE GAUCHE
				u1(i,Ny/3+5:Ny/3+15,1:24) = 0
				v1(i,Ny/3+5:Ny/3+15,1:24) = 0

				u1(i+1,Ny/3+5:Ny/3+15,1:24) = 0
                                v1(i+1,Ny/3+5:Ny/3+15,1:24) = 0

				u2(i,Ny/3+5:Ny/3+15,1:24) = 0
                                v2(i,Ny/3+5:Ny/3+15,1:24) = 0

                                u2(i+1,Ny/3+5:Ny/3+15,1:24) = 0
                                v2(i+1,Ny/3+5:Ny/3+15,1:24) = 0


				! COTE DROITE
				u1(i,Ny/3+5:Ny/3+15,26:50) = 0
				v1(i,Ny/3+5:Ny/3+15,26:50) = 0

				u1(i+1,Ny/3+5:Ny/3+15,26:50) = 0
                                v1(i+1,Ny/3+5:Ny/3+15,26:50) = 0

				u2(i,Ny/3+5:Ny/3+15,26:50) = 0
                                v2(i,Ny/3+5:Ny/3+15,26:50) = 0

                                u2(i+1,Ny/3+5:Ny/3+15,26:50) = 0
                                v2(i+1,Ny/3+5:Ny/3+15,26:50) = 0

			end if


		end do
	
	end subroutine


end program
