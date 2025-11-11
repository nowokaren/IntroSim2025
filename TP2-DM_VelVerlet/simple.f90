! Introducción a la Simulación Computacional
! Edición: 2025
! Docentes: Joaquín Torres y Claudio Pastorino

program simple 

    use ziggurat
    implicit none
    logical :: es
    integer :: seed, N, i, j, n_step, N_steps, N_minE, N_eq 
    real(8) :: L, dc, Vc, sigma, eps, dt, m, Vtot, Ktot, Etot, vpair
    real(8) :: d, dx, dy, dz, d2, d6, d12, fx, fy, fz, fmag
    real(8) :: rho, T_target, temp, scale
    real(8), allocatable :: r(:,:), v(:,:), f(:,:), E(:), T(:), sumv(:)
    character(80) :: filename


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

    N         = 108                    ! 108 partículas (FCC 3x3x3)
    sigma     = 1.0d0
    L         = ((N/0.4d0)**(1.0d0/3.0d0))*sigma  ! ρ = N/V = 0.4
    eps       = 1.0d0
    m         = 1.0d0
    dc        = 2.5d0*sigma
    T_target  = 1.5d0                  ! T = 1.5 ε/kB ; T = m <v²> / (3Nk_B)
    dt        = 0.005d0                ! paso temporal
    N_minE    = 10000                ! min energía
    N_eq      = 50000                    ! equilibración
    N_steps   = N_eq + N_minE
    
    allocate(r(N,3), v(N,3), f(N,3), E(N_steps), T(N_steps), sumv(3))

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
    temp = sum(v**2) / (3.0d0 * N)
    T = temp
    
    ! MINIMIZACIÓN DE ENERGÍA
    do n_step = 1, N_minE
        if (n_step > 1) then
            r = r + f*(dt*dt)/(2.0d0*m)   ! r(t+dt) = r(t) + f dt²/2m
            r = r - L * floor(r / L)
        end if
        Vtot = 0.0d0
        f = 0.0d0
        do i = 1, N-1
            do j = i+1, N
                dx = r(i,1)-r(j,1); dy = r(i,2)-r(j,2); dz = r(i,3)-r(j,3)
                dx = dx - L*nint(dx/L)
                dy = dy - L*nint(dy/L)
                dz = dz - L*nint(dz/L)
                d2 = dx*dx + dy*dy + dz*dz
                if (d2 >= dc**2 .or. d2 < 1d-12) cycle
                d = sqrt(d2)
                d6  = d2**3
                d12 = d6*d6
                vpair = 4.0d0*eps*(1.0d0/d12 - 1.0d0/d6) - Vc
                Vtot = Vtot + vpair
                fmag = 24.0d0*eps * (2.0d0*d12**(-1) - d6**(-1)) / d  
                fx = fmag*(dx/d); fy = fmag*(dy/d); fz = fmag*(dz/d)
                f(i,1) = f(i,1) + fx; f(j,1) = f(j,1) - fx
                f(i,2) = f(i,2) + fy; f(j,2) = f(j,2) - fy
                f(i,3) = f(i,3) + fz; f(j,3) = f(j,3) - fz
                !print*, i, j, d, vpair, fx, fy, fz
            end do
        end do
        E(n_step) = Vtot
        !print*, "step =", n_step, "Epot =", Vtot

        if (n_step == 1) then
            open(30, file='traj.vtf', status='replace')
            write(30,'(A)') '# LJ - IntroSim 2025'
            ! átomos
            do i = 1, N
                write(30,'(A,I0,A)') 'atom ',i-1,' name Ar'
                write(30,'(A,I0,A,F6.3)') 'atom ',i-1,' radius ',0.5
            end do
            ! c aja
            write(30,'(A,3F12.6)') 'unitcell ', L, L, L
        else
            open(30, file='traj.vtf', status='old', position='append')
        end if

        write(30,'(A)') 'timestep ordered'
        do i = 1, N
            write(30,'(3F16.8)') r(i,1), r(i,2), r(i,3)
        end do
        close(30)

    end do
    dt=dt/4
    ! TERMALIZACIÓN/equilibración
    do n_step = N_minE+1, N_steps
        v = v + f*dt/(2.0d0*m)  ! v(t+dt/2)
        r = r + v * dt          ! r(t+dt)
        r = r - L * floor(r / L)
        Vtot = 0.0d0; f = 0.0d0
        do i = 1, N-1
            do j = i+1, N
                dx = r(i,1)-r(j,1); dy = r(i,2)-r(j,2); dz = r(i,3)-r(j,3)
                dx = dx - L*nint(dx/L)
                dy = dy - L*nint(dy/L)
                dz = dz - L*nint(dz/L)
                d2 = dx*dx + dy*dy + dz*dz
                if (d2 >= dc**2 .or. d2 < 1d-12) cycle
                d = sqrt(d2)
                d6  = d2**3
                d12 = d6*d6
                vpair = 4.0d0*eps*(1.0d0/d12 - 1.0d0/d6) - Vc
                Vtot = Vtot + vpair
                fmag = 24.0d0*eps * (2.0d0*d12**(-1) - d6**(-1)) / d  
                fx = fmag*(dx/d); fy = fmag*(dy/d); fz = fmag*(dz/d)
                f(i,1) = f(i,1) + fx; f(j,1) = f(j,1) - fx
                f(i,2) = f(i,2) + fy; f(j,2) = f(j,2) - fy
                f(i,3) = f(i,3) + fz; f(j,3) = f(j,3) - fz
                !print*, i, j, d, vpair, fx, fy, fz
            end do
        end do
        E(n_step) = Vtot

        v = v+ f*dt/(2.0d0*m)   ! v(t+dt)
        temp = sum(v**2) /(3.0d0*N)
        T(n_step) = temp
        !print*, "step =", n_step, "Epot =", Vtot, " T =", temp

        if (n_step == 1) then
            open(30, file='traj.vtf', status='replace')
            write(30,'(A)') '# LJ - IntroSim 2025'
            ! átomos
            do i = 1, N
                write(30,'(A,I0,A)') 'atom ',i-1,' name Ar'
                write(30,'(A,I0,A,F6.3)') 'atom ',i-1,' radius ',0.5
            end do
            ! caja
            write(30,'(A,3F12.6)') 'unitcell ', L, L, L
        else
            open(30, file='traj.vtf', status='old', position='append')
        end if

        write(30,'(A)') 'timestep ordered'
        do i = 1, N
            write(30,'(3F16.8)') r(i,1), r(i,2), r(i,3)
        end do
        close(30)
    end do


    ! GUARDAR ENERGIA y Temperatura
    open(unit=21, file="equilibrio.dat", status='replace')
    write(21, *) '# i   E   T'
    do i = 1, N_steps
        write(21, '(I5, 2ES15.6E2)') i, E(i), T(i)
    end do
    close(21)

    

!! FIN FIN edicion
!! 
![No TOCAR]
! Escribir la última semilla para continuar con la cadena de numeros aleatorios 

        open(unit=10,file='seed.dat',status='unknown')
        seed = shr3() 
         write(10,*) seed
        close(10)
![FIN no Tocar]        


end program simple
