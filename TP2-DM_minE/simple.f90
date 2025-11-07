! Introducción a la Simulación Computacional
! Edición: 2025
! Docentes: Joaquín Torres y Claudio Pastorino

program simple 


    use ziggurat
    implicit none
    logical :: es
    integer :: seed,N, i, j,  cont, N_steps, n_step
    real (kind=8) :: d, d2, d6, d12, fx, fy, fz, dx, dy, dz, vpair
    real (kind=8) :: L, dc, Vc, sigma, eps, fmag, dt, m, Vtot, f_proj
    real (kind=8), allocatable  :: y(:), c(:,:), r(:,:),f(:,:), r_cont(:,:,:), Vpot(:), r_j(:), E(:)
    character(len=30) :: filename



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

    N_steps = 20000
    dt = 0.01d0 
    m = 1.0d0   
    N = 20
    L = 5.0d0
    eps=1.0d0
    sigma=1.0d0
    dc = 2.5*sigma
    
    allocate(r(N,3),f(N,3), r_cont(N,N,3),  E(N_steps))

    do i = 1, N
        do j = 1, 3
            r(i,j) = uni()*L
        end do
    end do

    do n_step = 1, N_steps
        print *, "step", n_step
        ! --- Condiciones periódicas ---
        if (n_step > 1) then
            r(:,:) = r(:,:) + f(:,:)*(dt*dt)/(2.0d0*m)   ! r(t+dt) = r(t) + f dt²/2m
            where(r >= L)  r = r - L
            where(r <  0)  r = r + L
        end if

        Vc = 4.0d0*eps*( (sigma/dc)**12 - (sigma/dc)**6 )
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

                ! energía (shifted)
                vpair = 4.0d0*eps*(1.0d0/d12 - 1.0d0/d6) - Vc
                Vtot = Vtot + vpair

                ! fuerza
                fmag = 24.0d0*eps * (2.0d0*d12**(-1) - d6**(-1)) / d  
                fx = fmag*(dx/d); fy = fmag*(dy/d); fz = fmag*(dz/d)

                ! 3ra ley
                f(i,1) = f(i,1) + fx; f(j,1) = f(j,1) - fx
                f(i,2) = f(i,2) + fy; f(j,2) = f(j,2) - fy
                f(i,3) = f(i,3) + fz; f(j,3) = f(j,3) - fz
                !print*, i, j, d, vpair, fx, fy, fz

                dx = r(1,1)-r(2,1); dy = r(1,2)-r(2,2); dz = r(1,3)-r(2,3)
                dx = dx - L*nint(dx/L); dy = dy - L*nint(dy/L); dz = dz - L*nint(dz/L)
                d = sqrt(dx*dx + dy*dy + dz*dz)
                if (d > 1d-10) then
                    f_proj = (f(1,1)*dx + f(1,2)*dy + f(1,3)*dz) / d   ! F proyectada
                    !print*, n_step, d, f_proj
                end if
            end do
        end do

        E(n_step) = Vtot
        print*, "Epot total =", Vtot, "  Media por partícula =", Vtot/N

        write(filename, '("estados/estado_", I0, ".dat")') n_step
        open(unit=20, file=filename, status='replace')
        write(20, *) '# i   rx ry rz   fx fy fz'
        do i = 1, N
            write(20, '(I5, 6ES15.6E2)') i, r(i,1), r(i,2), r(i,3), &
                                        f(i,1), f(i,2), f(i,3)
        end do
        close(20)  

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

    open(unit=21, file="energy.dat", status='replace')
    write(21, *) '# i   E'
    do i = 1, N_steps
        write(21, '(I5, 1ES15.6E2)') i, E(i)
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
