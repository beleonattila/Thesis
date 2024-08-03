#!/usr/bin/bash

dir="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/alphafold" # change this folder name!!!!
logdir="$dir/logs"
datadir="$dir/data"
resultdir="$dir/results"

for file in ${datadir}/*.faa ; do
  if [ ! -f $file ] ; then
          echo "File does not exits:" $f
          exit
  fi

        # save file name
#  filename=`echo $file | cut -d '/' -f 9 | cut -d '.' -f 1`
  filename=`echo $file | cut -d '/' -f 10 | cut -d '.' -f 1`
  
        #echo $file
        echo $filename

  # Make script for each sample
  script="${logdir}/${filename}_run.sh"
        echo $script
  # Specify bash shebang
  printf "#!/usr/bin/bash\n" > $script 

  # Introduction to script
  printf "\n# Script generated the: $(date +'%Y-%m-%d')\n" >> $script

  # Load required modules
  printf "\n# Load required modules\n" >> $script
  printf "module load tools cuda/toolkit/11.4.1 cudnn/11.4-8.2.4.15 tensorrt/11.4-8.2.1.8 anaconda3/2023.03 perl/5.36.1 gcc/10.3.0 openmpi/4.1.6 hhsuite/3.3.0 hmmer/3.3.2 kalign/3.3.1 cuda/toolkit/11.4.1 cudnn/11.4-8.2.4.15 tensorrt/11.4-8.2.1.8 openmpi/4.1.6 hhsuite/3.3.0 hmmer/3.3.2 kalign/3.3.1 alphafold/2.3.1\n" >> $script
 
  # Run AlphaFold
  printf "\n# Run AlphaFold\n" >> $script
  printf "/usr/bin/time -v -o ${logdir}/${filename}.time /services/tools/alphafold/2.3.2/run_alphafold.sh -d /home/databases/alphafold -o $resultdir -f $file -t 2024-01-01 -l 1 -r false\n\n" >> $script

  # Qsub script
  qsub -W group_list=cge -A cge -d `pwd` -l nodes=1:gpus=1:ppn=40,mem=150gb,walltime="10:00:00" -r y -N $filename -o $logdir -e $logdir $script
  sleep 0.5

done
