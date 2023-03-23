#! /bin/bash
# ============================================================
# bash to merge bam files for each individual with more than one fastq.gz files, work with parallel
# ============================================================

# get command line options
# ============================================================

args=`getopt -o "r:,l:,t:,o:,i:,b:" -l "resource:,list:,tmp:,out:,individual:,bamInd:" -- "$@"`
>&2 echo "command arguments given: $args"
eval set -- "$args"

tmp=/tmp/

# parse arguments
# ============================================================

while true;
do
  case $1 in

    -r|--resource)
      resource=$2
      shift 2;;

    -l|--list)
			list=$2
			shift 2;; 

  	-t|--tmp)
			tmp=$2
			shift 2;;

    -o|--out)
			out=$2
			shift 2;;

    -i|--individual)
		  individual=$2
			shift 2;;
    
    -b|--bamInd)
		  bamInd=$2
			shift 2;;

    --)
      shift
      break;;

  esac
done

# prepare directory
# ============================================================

if [[ -e $out/$individual ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: $out/$individual exists."
  exit 1
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: makedir $out/$individual."
  mkdir $out/$individual
	mkdir $out/$individual/log
fi

# check for file availability
# ============================================================

if [[ -e $resource/resource.env && -e $resource/sbatch/mergeBamInd.sbatch ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: all slurm scripts found, starting to merge"
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: cannot find some slurm scripts."
  exit 1
fi

# merge sort bam per individual
# ============================================================

jobID=$(sbatch --export=env=$resource/resource.env,dir=$out/$individual,tmp=$tmp,list=$list,bamInd=$bamInd --output=$out/$individual/merge.out --error=$out/$individual/merge.err $resource/sbatch/mergeBamInd.sbatch | cut -d " " -f 4)
sleep 5m
jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
while [ "$jobState" == "PD" ] || [ "$jobState" == "R" ] || [ "$jobState" == "CG" ]
do
	sleep 5m
	jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
done

# check if run success
reRunCounter=0

until [ `wc -l $out/$individual/merge.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" $out/$individual/merge.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 2 ]
do
	rm -r $out/$individual/*
	mkdir $out/$individual/log
	#sacct --units=G --format="User,state%15,Account%15,JobID%20,Jobindividual%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqGRES,AllocGRES" -j $jobID >> $out/$individual/merge.fail.out
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting merge job."
	mergeTime=$(($reRunCounter * 4 + 12))
	jobID=$(sbatch --export=env=$resource/resource.env,dir=$out/$individual,tmp=$tmp,list=$list,bamInd=$bamInd --time=$mergeTime:00:00 --output=$out/$individual/merge.out --error=$out/$individual/merge.err $resource/sbatch/mergeBamInd.sbatch | cut -d " " -f 4)
	sleep 5m
	jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
	while [ "$jobState" == "PD" ] || [ "$jobState" == "R" ] || [ "$jobState" == "CG" ]
	do
		sleep 5m
		jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l $out/$individual/merge.err | awk '{print $1}'` == 0 ] && [ `grep "done.main.process" $out/$individual/merge.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed merge alignments."
	rm -r `echo $bamInd | tr "," "\n" | grep -w -f - $list | awk '{print "'$tmp'/"$1""}'`
	#sacct --units=G --format="User,state%15,Account%15,JobID%20,Jobindividual%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqGRES,AllocGRES" -j $jobID >> $out/$individual/merge.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):error: error for merge."
	exit 1
fi
