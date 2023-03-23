#! /bin/bash
# ============================================================
# bash to call sbatch hisat2 build and hisat2 mapping and work with parallel
# ============================================================

# get command line options
# ============================================================

args=`getopt -o "r:,w:,o:,i:,l:,b:" -l "resource:,wdir:,out:,individual:,line:,bam:" -- "$@"`
>&2 echo "command arguments given: $args"
eval set -- "$args"

# parse arguments
# ============================================================

while true;
do
  case $1 in

    -r|--resource)
      resource=$2
      shift 2;;

    -w|--wdir)
			wdir=$2
			shift 2;; 

    -o|--out)
			out=$2
			shift 2;;
			
    -i|--individual)
		  individual=$2
			shift 2;;

    -l|--line)
		  line=$2
			shift 2;;
    
    -b|--bam)
      bam=$2
      shift 2;;
			
    --)
      shift
      break;;

  esac
done

# check for file availability
# ============================================================

if [[ -e $resource/resource.env && -e $resource/sbatch/stringtie.sbatch ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: all slurm scripts found."
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: cannot find some slurm scripts."
  exit 1
fi

if [[ -e $bam ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info:starting stringtie"
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):info:cannot find $bam file"
  exit 1
fi

# stringtie 
# ============================================================

mkdir "$out"/"$individual"

jobID=$(sbatch --export=env=$resource/resource.env,dir=$out/$individual,gtf=$wdir/$line/lift.gtf,bam=$bam --output="$out"/"$individual"/stringtie.out --error="$out"/"$individual"/stringtie.err $resource/sbatch/stringtie.sbatch | cut -d " " -f 4)

sleep 1m
jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
while [ "$jobState" == "PD" ] || [ "$jobState" == "R" ] || [ "$jobState" == "CG" ]
do
  sleep 5m
  jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
done

reRunCounter=0
until [ `grep -i "error" "$out"/"$individual"/stringtie.err | wc -l | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" "$out"/"$individual"/stringtie.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 2 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting stringtie job."
	#sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqGRES,AllocGRES" -j $jobID >> $out/$individual/stringtie.fail.out
	buildTime=$(($reRunCounter * 4 + 8))
	buildMem=$(($reRunCounter * 4 + 20))G

  jobID=$(sbatch --export=env=$resource/resource.env,dir=$out/$individual,gtf=$wdir/$line/lift.gtf,bam=$bam --time=$buildTime:00:00 --mem=$buildMem --output="$out"/"$individual"/stringtie.out --error="$out"/"$individual"/stringtie.err $resource/sbatch/stringtie.sbatch | cut -d " " -f 4)  
	sleep 5m
	jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
	while [ "$jobState" == "PD" ] || [ "$jobState" == "R" ] || [ "$jobState" == "CG" ]
	do
		sleep 5m
		jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `grep -i "error" "$out"/"$individual"/stringtie.err | wc -l | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" "$out"/"$individual"/stringtie.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed stringtie."
	#sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqGRES,AllocGRES" -j $jobID >> "$out"/"$individual"/stringtie.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):info: error for stringtie."
	exit 1
fi
