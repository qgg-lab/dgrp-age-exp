#! /bin/bash
# ============================================================
# bash to call sbatch hisat2 mapping and merge, work with parallel
# ============================================================

# get command line options
# ============================================================

args=`getopt -o "r:,w:,t:,o:,l:,f:,s:,i:,c:,a:" -l "resource:,wdir:,tmp:,out:,line:,fastq:,sample:,individual:,fc:,lane:" -- "$@"`
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

    -w|--wdir)
			wdir=$2
			shift 2;;

  	-t|--tmp)
			tmp=$2
			shift 2;;

    -o|--out)
			out=$2
			shift 2;;

    -l|--line)
		  line=$2
			shift 2;;
			
    -f|--fastq)
		  fastq=$2
			shift 2;;
			
    -s|--sample)
		  sample=$2
			shift 2;;

    -i|--individual)
		  individual=$2
			shift 2;;
	
    -c|--fc)
		  fc=$2
			shift 2;;

    -a|--lane)
		  lane=$2
			shift 2;;
			
    --)
      shift
      break;;

  esac
done

# prepare directory
# ============================================================

if [[ -e $tmp/$sample ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: $tmp/$sample exists."
  exit 1
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: makedir $tmp/$sample."
  mkdir $tmp/$sample
fi

if [[ -e $out/$sample ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: $out/$sample exists."
  exit 1
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: makedir $out/$sample."
  mkdir $out/$sample
	mkdir $out/$sample/log
fi

# check for file availability
# ============================================================

if [[ -e $resource/resource.env && -e $resource/sbatch/hisat2.mapping.sbatch ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: all slurm scripts found."
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: cannot find some slurm scripts."
  exit 1
fi


if [ `ls $wdir/$line | grep "ht2" | wc -l | awk '{print $1}'` -eq 8 ] && [ `find $wdir/$line -type f -size 0b | grep "ht2" | wc -l | awk '{print $1}'` -eq 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: starting hisat2 mapping"
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: cannot find index"
  exit 1
fi

# =====================
# = 1. hisat2.mapping =
# =====================

jobID=$(sbatch --export=env="$resource"/resource.env,fastq="$fastq",dir="$out"/"$sample",individual="$individual",fc="$fc",lane="$lane",index="$wdir"/"$line"/lift.genome,tmp="$tmp"/"$sample"/ --output="$out"/"$sample"/hisat2.mapping.out --error="$out"/"$sample"/hisat2.mapping.err "$resource"/sbatch/hisat2.mapping.sbatch | cut -d " " -f 4)
sleep 5m
jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
while [ "$jobState" == "PD" ] || [ "$jobState" == "R" ] || [ "$jobState" == "CG" ]
do
	sleep 5m
	jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
done

# check if run success
reRunCounter=0

until [ `wc -l "$out"/"$sample"/hisat2.mapping.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" "$out"/"$sample"/hisat2.mapping.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 2 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting hisat2.mapping job."
	rm $tmp/$sample/*
	#sacct --units=G --format="User,state%15,Account%15,jobID%20,jobID%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqGRES,AllocGRES" -j $jobID >> $out/$sample/hisat2.mapping.fail.out
	mappingTime=$(($reRunCounter * 4 + 12))
	jobID=$(sbatch --export=env="$resource"/resource.env,fastq="$fastq",dir="$out"/"$sample",individual="$individual",fc="$fc",lane="$lane",index="$wdir"/"$line"/lift.genome,tmp="$tmp"/"$sample"/ --time=$mappingTime:00:00 --output="$out"/"$sample"/hisat2.mapping.out --error="$out"/"$sample"/hisat2.mapping.err "$resource"/sbatch/hisat2.mapping.sbatch | cut -d " " -f 4)
	sleep 5m
	jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
	while [ "$jobState" == "PD" ] || [ "$jobState" == "R" ] || [ "$jobState" == "CG" ]
	do
		sleep 5m
		jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l "$out"/"$sample"/hisat2.mapping.err | awk '{print $1}'` == 0 ] && [ `grep "done.main.process" "$out"/"$sample"/hisat2.mapping.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed hisat2.mapping."
  #sacct --units=G --format="User,state%15,Account%15,jobID%20,jobID%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqGRES,AllocGRES" -j $jobID >> $out/$sample/hisat2.mapping.out
else
	>&2 echo "$(date):error: error for hisat2.mapping."
	exit 1
fi
