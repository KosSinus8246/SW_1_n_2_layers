program SW1L
	implicit none
	
	real(8) :: pi, omega, g, Lx, Ly, dx, dy, dt, time, f
	integer :: i, j, k, Nx, Ny, Nt, tspawn

	character(50) :: IC, obstacle
	logical :: rotating_frame 

	real(8), allocatable :: x(:), y(:), t(:), u(:,:,:), v(:,:,:), eta(:,:,:), H(:,:)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	pi = 3.14159265359
	omega = 2*pi/(24*3600)	
	g = 9.81
	tspawn = 2

	IC = 'eta_detroit'
	obstacle = 'ile'
	rotating_frame = .false.


	Lx = 5*1e3
	Ly = 5*1e3

	Nx = 50
	Ny = 50

	dx = Lx/Nx
	dy = Ly/Ny

	dt = 20
	time = 5000
	Nt = time/dt

	allocate(t(Nt), x(Nx), y(Ny))

	do i=1,Nt-1
		t(i+1) = t(i) + dt
	end do
	do j=1,Nx-1
		x(j+1) = x(j) + dx
	end do
	do k=1,Ny-1
                y(k+1) = y(k) + dy
        end do

	allocate(u(Nt, Nx, Ny), v(Nt, Nx, Ny), eta(Nt, Nx, Ny), H(Nx, Ny))


	if (rotating_frame .eqv. .false.) then
		f = 0.
		print *, 'no coriolis'
		call integrator_LP(Nt, Nx, Ny, u, v, eta, dx, dy, dt, g, H, f, obstacle, tspawn)
	elseif (rotating_frame .eqv. .true.) then
		print *, 'bz ta mer'
		f = 2*omega*sin(45.)
		call integrator_LP(Nt, Nx, Ny, u, v, eta, dx, dy, dt, g, H, f, obstacle, tspawn)

	end if


	print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	print *, 'OK'


contains

	subroutine integrator_LP(Nt, Nx, Ny, u, v, eta, dx, dy, dt, g, H, f, obstacle, tspawn)
		
		!integer :: i, j, k

		integer, intent(in) :: Nt, Nx, Ny, tspawn
		real(8), intent(in) :: dx, dy, dt, g, f
		
		real(8), intent(out) :: u(:,:,:), v(:,:,:), eta(:,:,:), H(:,:)
		
		character(50), intent(in) :: obstacle

		do i=2,Nt-1
			do j=2,Nx-1
				do k=2,Ny-1
					u(i,1,:) = 0.
					u(i,Nx,:) = 0.
					v(i,:,1) = 0.
					v(i,:,Ny) = 0.
		
					if (obstacle == 'detroit') then
						! COTE GAUCHE
						u(i,Ny/3+5:Ny/3+15,1:24) = 0.
						v(i,Ny/3+5:Ny/3+15,1:24) = 0.

						! COTE DROITE
						u(i,Ny/3+5:Ny/3+15,26:50) = 0.
						v(i,Ny/3+5:Ny/3+15,26:50) = 0.

					elseif (obstacle == 'ile') then
						u(i,Ny/3+5:Ny/3+15,20:30) = 0.
						v(i,Ny/3+5:Ny/3+15,20:30) = 0.
					end if

					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					!SCHEMA MIXTE : 1) Euler ; 2) LP ; 3) Euler ; 4) LP.
					
					if (t(i) <= tspawn) then
						!print('Euler init')
						u(i+1,j,k) = u(i,j,k)-dt*(g* (eta(i,j+1,k)-eta(i,j-1,k))/(2*dx) + f*v(i,j,k))
						v(i+1,j,k) = v(i,j,k)-dt*(g* (eta(i,j,k+1)-eta(i,j,k-1))/(2*dy) - f*u(i,j,k))
						eta(i+1,j,k) = eta(i,j,k)-dt*(H(j,k)*((u(i,j+1,k)-u(i,j-1,k))/(2*dx) + (v(i,j,k+1)-v(i,j,k-1))/(2*dy)))

					elseif (t(i) > tspawn) then 
						!print('Leap-Frog')
						u(i+1,j,k) = u(i-1,j,k)-2*dt*(g* (eta(i,j+1,k)-eta(i,j-1,k))/(2*dx) + f*v(i,j,k))
						v(i+1,j,k) = v(i-1,j,k)-2*dt*(g* (eta(i,j,k+1)-eta(i,j,k-1))/(2*dy) - f*u(i,j,k))
						eta(i+1,j,k) = eta(i-1,j,k)-2*dt*(H(j,k)*((u(i,j+1,k)-u(i,j-1,k))/(2*dx) + (v(i,j,k+1)-v(i,j,k-1))/(2*dy)))
					end if


				end do
			end do
		end do
	
	end subroutine


end program
