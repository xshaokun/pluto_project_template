#BSUB -J p2r003-thin-coarse-trc
#BSUB -q astron
#BSUB -n 100
#BSUB -o ./p2r003-trc/p2r003-thin-coarse-trc.out
#BSUB -e ./p2r003-trc/p2r003-thin-coarse-trc.err
#BSUB -R span[ptile=25]
mpirun -np 100 ./pluto -restart 12 -x3jet
