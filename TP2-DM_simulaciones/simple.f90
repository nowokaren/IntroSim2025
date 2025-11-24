
    ! Introducción a la Simulación Computacional
    ! Edición: 2025
    ! Docentes: Joaquín Torres y Claudio Pastorino

    program simple 

        use ziggurat
        implicit none
        logical :: es
        integer :: seed, N, i, j, n_step, N_steps, N_minE, N_eq, ios, nargs
        real(8) :: L, dc, Vc, sigma, eps, dt, m
        real(8) :: gamma, kB, T_target, potential, virial, rho
        real(8) :: kinetic, total_energy, temp, pressure, lgv_term
        real(8), allocatable :: r(:,:), v(:,:), f(:,:), sumv(:)
        !real(8), allocatable :: K(:), Vpot(:), E(:), T(:), P(:), Vir(:)   ! K, V, E, T, P = ARRAYS
        character(80) :: filename, dirname
        character(len=20) :: arg1, arg2

        integer :: dump_freq = 500    ! Frecuencia para dump traj y mediciones

        real(8) ::  sum_p = 0.d0, sum_p2 = 0.d0, avg_p, sigma_p  ! Para promedios on-fly
        integer :: count_meas = 0


        T_target = 1.1d0
        rho      = 0.4d0
        arg1 = ''
        arg2 = ''
        call get_command_argument(1, arg1)
        call get_command_argument(2, arg2)
        if (len_trim(arg1) > 0) read(arg1, *) T_target
        if (len_trim(arg2) > 0) read(arg2, *) rho


    ![NO TOCAR] Inicializa generador de número random

        inquire(file='seed.dat',exist=es)
        if(es) then
            open(unit=10,file='seed.dat',status='old')
            read(10,*) seed
            close(10)
            !print *,"  * Leyendo semilla de archivo seed.dat"
        else
            seed = 24583490
        end if

        call zigset(seed)
    ![FIN NO TOCAR]    

        N         = 200                   
        sigma     = 1.0d0
       !rho       = arg2
        L         = ((N/rho)**(1.0d0/3.0d0))*sigma  ! ρ = N/V = 0.4
        eps       = 1.0d0
        m         = 1.0d0
        kB        = 1.0d0
        dc        = 2.5d0*sigma
       !T_target  = arg1                  ! T = 1.5 ε/kB ; T = m <v²> / (3Nk_B)
        gamma     = 0.3d0                  ! Coeficiente de fricción
        dt        = 0.002d0               ! paso temporal
        N_minE    = 30000                     ! min energía
        N_eq      = 2000000                ! equilibración
        N_steps   = N_eq + N_minE
        lgv_term = sqrt(2.0d0 * gamma * m * kB * T_target / dt)

  
        write(dirname, '(A, I0, A, F4.2, A, F5.3, A, F5.3)') &
                'simulaciones/N', N, '_T', T_target, '_rho', real(N / (L**3)), '_dt', dt

        call system('mkdir -p ' // trim(dirname))

        
        allocate(r(N,3), v(N,3), f(N,3), sumv(3))
        !allocate(K(N_steps), Vpot(N_steps), E(N_steps), T(N_steps), P(N_steps), Vir(N_steps))
        Vc = 4.0d0*eps*( (sigma/dc)**12 - (sigma/dc)**6 )

        ! POSICIONES
        do i = 1, N
            do j = 1, 3
                r(i,j) = uni()*L
            end do
        end do

        ! VELOCIDADES
        do i = 1, N
            do j = 1, 3
                v(i,j) = rnor()
            end do
        end do

        sumv = sum(v, dim=1)
        v = v - spread(sumv/N, 1, N) 
        temp = sum(v**2) / (3.0d0 * N)  
        v = v * sqrt(T_target / temp)
        call force(r, f, potential, virial, L, N, dc, eps, sigma, Vc)

        dt = 0.0001d0
                
        ! MINIMIZACIÓN DE ENERGÍA (NVE)
        do n_step = 1, N_minE
            if (n_step > 1) then
                r = r + f*(dt*dt)/(2.0d0*m)   ! r(t+dt) = r(t) + f dt²/2m
                r = r - L * floor(r / L)
            end if
            call force(r, f, potential, virial, L, N, dc, eps, sigma, Vc)
        end do
        
        ! Guardar energia, temperatura y presion
        open(21, file=trim(dirname)//'/mediciones.dat', status='replace')
        write(21,*) '# step     K          V          E         T         P         Vir'
        
        open(30, file=trim(dirname)//'/traj.vtf', status='replace')
        do i = 1, N
            write(30,'(A,I0,A,F6.3)') 'atom ',i-1,' name Ar radius ',0.5
        end do
        write(30,'(A,3F12.6)') 'unitcell ', L, L, L

        dt = 0.002d0    
        ! TERMALIZACIÓN/equilibración
        do n_step = N_minE+1, N_steps

            call verlet_positions(r, v, f, dt, m, L, N)
            call force(r, f, potential, virial, L, N, dc, eps, sigma, Vc)
            call lgv_force(v, f, N, m, dt, T_target, gamma, lgv_term)
            call verlet_velocities(v, f, dt, m, N)
            call measurements(r, v, potential, virial, &
                            kinetic, total_energy, temp, pressure, &
                            N, m, L)

            count_meas = count_meas + 1
            sum_p = sum_p + pressure
            sum_p2 = sum_p2 + pressure**2

            if (mod(n_step, dump_freq) == 0) then
                write(21,'(I10,6ES14.6)') n_step, kinetic, potential, total_energy, temp, pressure, virial
                write(30,'(A)') 'timestep ordered'
                do i = 1, N
                    write(30,'(3E16.8)') r(i,1), r(i,2), r(i,3)
                end do
            end if

        end do

        close(21)
        close(30)
        
        if (count_meas > 0) then
            avg_p = sum_p / real(count_meas,8)
            sigma_p = sqrt( sum_p2 / real(count_meas,8) - avg_p**2 )
            open(23, file=trim(dirname)//'/summary.dat', status='replace')
            write(23,*) '# Avg_P     Sigma_P'
            write(23,'(2ES14.6)') avg_p, sigma_p
            close(23)
        end if

        
        ! Guardar configuración 
        open(unit=99, file=trim(dirname)//'/config.txt', status="replace")
        write(99,'(A, I0)') "N = ", N
        write(99,'(A, F12.6)') "L = ", L
        write(99,'(A, F12.6)') "sigma = ", sigma
        write(99,'(A, F12.6)') "eps = ", eps
        write(99,'(A, F12.6)') "m = ", m
        write(99,'(A, F12.6)') "dc = ", dc
        write(99,'(A, F12.6)') "dt = ", dt
        write(99,'(A, F12.6)') "gamma = ", gamma
        write(99,'(A, F12.6)') "kB = ", kB
        write(99,'(A, F12.6)') "T_target = ", T_target
        write(99,'(A, I0)') "N_minE = ", N_minE
        write(99,'(A, I0)') "N_eq = ", N_eq
        write(99,'(A, I0)') "N_steps = ", N_steps
        write(99,'(A, F12.6)') "Vc = ", Vc
        write(99,'(A, F12.6)') "rho = ", N / (L**3)
        call system('date "+Fecha: %Y-%m-%d %H:%M:%S" >> ' // trim(dirname)//'/config.txt')

        close(99)


    !! FIN FIN edicion
    !! 
    ![No TOCAR]
    ! Escribir la última semilla para continuar con la cadena de numeros aleatorios 

            open(unit=10,file='seed.dat',status='unknown')
            seed = shr3() 
            write(10,*) seed
            close(10)
    ![FIN no Tocar]     
    
    contains

        function pbc(ri, rj, L) result(diff)
            real(8), intent(in) :: ri(3), rj(3), L
            real(8) :: diff(3), dr(3)
            integer :: k

            dr = ri - rj
            do k = 1, 3
                diff(k) = dr(k) - L * nint(dr(k) / L)
            end do
        end function pbc

    ! Calcula fuerzas LJ, energía potencial y virial con PBC
    subroutine force(r, f, potential, virial, L, N, dc, eps, sigma, Vc)
        real(8), intent(in) :: r(N,3), L, dc, eps, sigma, Vc
        integer, intent(in) :: N
        real(8), intent(out) :: f(N,3), potential, virial
        real(8) :: diff(3), d2, d, d6, d12, vpair, fmag
        integer :: i, j

        potential = 0.0d0
        virial = 0.0d0
        f = 0.0d0

        do i = 1, N-1
            do j = i+1, N
                diff = pbc(r(i,:), r(j,:), L)
                d2 = sum(diff**2)
                if (d2 >= dc**2 .or. d2 < 1d-12) cycle

                d = sqrt(d2)
                d6 = d2**3
                d12 = d6 * d6

                vpair = 4.0d0*eps*(1.0d0/d12 - 1.0d0/d6) - Vc
                potential = potential + vpair

                fmag = 24.0d0*eps * (2.0d0/d12 - 1.0d0/d6) / d
                f(i,:) = f(i,:) + fmag * (diff / d)
                f(j,:) = f(j,:) - fmag * (diff / d)

                virial = virial + dot_product(diff, fmag * (diff / d)) ! suma r_ij · f_ij
            end do
        end do
    end subroutine force

    ! Actualiza v(t + dt/2), r(t + dt) y aplica PBC
    subroutine verlet_positions(r, v, f, dt, m, L, N)
        integer, intent(in) :: N
        real(8), intent(in) :: dt, m, L
        real(8), intent(inout) :: r(N,3), v(N,3), f(N,3)
        integer :: i

        v = v + f * (dt / (2.0d0 * m))        !  v(t + dt/2)
        r = r + v * dt                        !  r(t + dt)        
        r = r - L * floor(r / L)              !  PBC

    end subroutine verlet_positions


    ! Actualiza velocidades finales v(t + dt)
    subroutine verlet_velocities(v, f, dt, m, N)
        integer, intent(in) :: N
        real(8), intent(in) :: dt, m
        real(8), intent(inout) :: v(N,3), f(N,3)

        v = v + f * (dt / (2.0d0 * m))
    end subroutine verlet_velocities


    ! Termostato de Langevin
    subroutine lgv_force(v, f, N, m, dt, T_target, gamma, lgv_term)
        integer, intent(in) :: N
        real(8), intent(in) :: m, dt, T_target, gamma, lgv_term
        real(8), intent(inout) ::f(N,3), v(N,3)
        integer :: i
        do i = 1, N
            ! Fricción
            f(i,:) = f(i,:) - gamma * m * v(i,:)
            
            ! Ruido
            f(i,1) = f(i,1) + lgv_term * rnor()
            f(i,2) = f(i,2) + lgv_term * rnor()
            f(i,3) = f(i,3) + lgv_term * rnor()
        end do
    end subroutine lgv_force

    ! Calcula cantidades físicas instantáneas
    subroutine measurements(r, v, potential, virial, kinetic, total_energy, &
                            temp, pressure, N, m, L)
        integer, intent(in) :: N
        real(8), intent(in) :: r(N,3), v(N,3), potential, virial, m, L
        real(8), intent(out) :: kinetic, total_energy, temp, pressure

        real(8) :: volume, density

        
        volume = L * L * L                    ! Volumen
        density = real(N, 8) / volume         ! Densidad
        kinetic = 0.5d0 * m * sum(v**2)       ! Energía cinética
        temp = sum(v**2) / (3.0d0 * N)        ! Temperatura instantánea
        total_energy = kinetic + potential    ! Energía total

        pressure = density * temp + virial / (3.0d0 * volume)

    end subroutine measurements


    end program simple