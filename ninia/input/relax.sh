f'''

#PBS -N {self.prefix}
#PBS -l mem={self.memory}gb,nodes=1:ppn={self.cpus},walltime={self.hours}:00:00
#PBS -q general

start=`date +%s`

mpirun -np {self.cpus} pw.x -nk 2 < {self.input_dir}/{self.prefix}.i > {self.input_dir}/{self.prefix}.out

end=`date +%s`

runtime=$((end-start))
echo $runtime
echo "Finished"

'''