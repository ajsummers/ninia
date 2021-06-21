f'''

#PBS -N Ru_r8_b1
#PBS -l mem=50gb,nodes=1:ppn=8,walltime=30:00:00
#PBS -q general

start=`date +%s`

mpirun -np 8 pw.x -nk 2 < /gpfs01/home/ajs0201/workQE/input/convergence-tests/round8/Ru_r8_b1.i > /gpfs01/home/ajs0201/workQE/input/convergence-tests/round8/Ru_r8_b1.out

end=`date +%s`

runtime=$((end-start))
echo $runtime
echo "Finished"

'''