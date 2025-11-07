! Introducción a la Simulación Computacional
! Edición: 2025
! Docentes: Joaquín Torres y Claudio Pastorino

program simple 


    use ziggurat
    implicit none
    logical :: es
    integer :: seed,N, i, j,  cont
    real (kind=8) :: aux, Vpot_tot, d, fx, fy, fz, dx, dy, dz, d_old, L, d_c, Vc, sigma, eps, fmag
    real (kind=8), allocatable  :: y(:), c(:,:), r(:,:),f(:,:), r_cont(:,:,:), Vpot(:), r_j(:)


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


    N = 10
    L = 1.0d0
    eps=1.0d0
    sigma=1.0d0
    d_c = 2.5*sigma

    
    allocate(r(N,3),f(N,3), r_cont(N,N,3), Vpot(N), r_j(3))
    r_cont=0.0
    do i = 1, N
        do j = 1, 3
            aux=uni()*L
            r(i,j) = aux
        end do
    end do

    Vc=4*eps*(-(sigma/d_c)**6 + (sigma/d_c)**12)
    do i = 1, N
        Vpot(i)=0.0
        do j = i+1, N
            !d = sqrt((r(i,1)-r(j,1))**2 + (r(i,2)-r(j,2))**2 + (r(i,3)-r(j,3))**2)
            
            dx = r(i,1) - r(j,1)
            dy = r(i,2) - r(j,2)
            dz = r(i,3) - r(j,3)
            dx = dx - L * nint(dx / L)  ! nint = nearest integer
            dy = dy - L * nint(dy / L)
            dz = dz - L * nint(dz / L)
            d = sqrt(dx**2 + dy**2 + dz**2)
            
            if (d > 1.0d-6) then
                Vpot(i) = Vpot(i) + 4*eps*((sigma/d)**12 - (sigma/d)**6) - Vc
            end if

        end do
    end do
    
    print*, "Potencial total:", sum(Vpot)


f = 0.0d0
do i = 1, N-1
    do j = i+1, N
        ! PBC 
        dx = r(i,1) - r(j,1)
        dy = r(i,2) - r(j,2)
        dz = r(i,3) - r(j,3)
        dx = dx - L * nint(dx / L)
        dy = dy - L * nint(dy / L)
        dz = dz - L * nint(dz / L)
        d = sqrt(dx*dx + dy*dy + dz*dz)

        !  CUTOFF  (particulas muy lejos)
        if (d > 1.0d-6) then
            !  F = -dV/dr * (dx/d, ...) 
            fmag = 24.0d0 * eps * (2.0d0*(sigma**12 / d**13) - (sigma**6 / d**7))
            fx = fmag * (dx / d)
            fy = fmag * (dy / d)
            fz = fmag * (dz / d)

            ! suma en i y resta en j 
            f(i,1) = f(i,1) + fx
            f(i,2) = f(i,2) + fy
            f(i,3) = f(i,3) + fz
            f(j,1) = f(j,1) - fx
            f(j,2) = f(j,2) - fy
            f(j,3) = f(j,3) - fz
        print*, i, j, fx, fy, fz  !print para comparar con el grafico: componente de fuerza de cada par
        end if
    end do
end do


open(unit=20, file='forces.dat', status='replace')
write(20, *) '# i   rx ry rz   fx fy fz   V   r_contx r_conty r_contz'
do i = 1, N
    write(20, '(I5, 10ES15.6E2)') i, r(i,1), r(i,2), r(i,3), &
                                  f(i,1), f(i,2), f(i,3), Vpot(i), &
                                  r_cont(i,1,1), r_cont(i,1,2), r_cont(i,1,3)
end do
close(20)







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
