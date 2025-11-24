#!/bin/bash
# Script para correr todas las simulaciones de la Guía 5

EXEC=./simple
mkdir -p simulaciones

#  1. Variando densidad (T = 1.1) 
T=1.1
echo "     Simulaciones a T = ${T} variando densidad "


#for rho in 0.001 0.01 0.1 0.8 0.9 1.0; do
for rho in 1.0 0.9 0.8 0.78 0.71 0.64 0.57 0.5  ; do
    echo "-> Ejecutando T=${T}, rho=${rho}"
    ${EXEC} ${T} ${rho}
done


# 10 puntos entre 0.1 y 0.8
!for rho in $(seq 0.15 0.07 0.8); do
!    printf -v rho_fmt "%.3f" "$rho"
!    rho="${rho/,/.}"
!    echo "-> Ejecutando T=${T}, rho=${rho}"
!    ${EXEC} ${T} ${rho}
!done


# 2. Variando temperatura (rho = 0.4)
!rho=0.4
!echo "---- Simulaciones a rho = ${rho} variando temperatura "

!for T in $(LC_NUMERIC=C seq 0.7 0.07 1.4); do
!    echo "-> Ejecutando T=${T}, rho=${rho}"
!    ${EXEC} ${T} ${rho}
!done

# 3. Distribución g(r) 
!echo "----   Simulaciones específicas para g(r) "

!declare -a Tvals=(1.1 1.4) 
!declare -a Rvals=(0.3 0.3)

!for i in ${!Tvals[@]}; do
!    T=${Tvals[$i]}
!    rho=${Rvals[$i]}
!    echo "-> Ejecutando T=${T}, rho=${rho}"
!    ${EXEC} ${T} ${rho}
!done
