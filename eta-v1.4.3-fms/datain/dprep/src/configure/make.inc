FC = ftn
CC = 
FFLAGS=  -static-nvidia -mcmodel=medium -Mlarge_arrays -DLITTLE -Mpreprocess
CFLAGS= -mcmodel=medium -Mlarge_arrays
CPP= pgcc -mcmodel=medium -Mlarge_arrays
CFLAGS = 