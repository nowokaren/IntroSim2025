! simple.f90 - Guía 1: Básicos de Fortran
program simple
  use ziggurat
  implicit none
  logical :: es
  integer :: seed, i, N
  real(8) :: x, y

  ! Leer N de input.dat
  open(10, file='input.dat', status='old')
  read(10,*) N
  close(10)

  ! Inicializar generador
  inquire(file='seed.dat', exist=es)
  if (es) then
     open(10, file='seed.dat', status='old')
     read(10,*) seed
     close(10)
  else
     seed = 123456789
  end if
  call zigset(seed)

  ! Generar N pares (x,y) y guardar si x < y
  open(20, file='output.dat', status='replace')
  do i = 1, N
     x = uni()
     y = uni()
     if (x < y) then
        write(20,*) x, y
     end if
  end do
  close(20)

  ! Guardar seed
  open(10, file='seed.dat', status='replace')
  seed = shr3()
  write(10,*) seed
  close(10)

end program simple
