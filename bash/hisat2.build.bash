#! /bin/bash
# ============================================================
# bash to call sbatch hisat2 build and hisat2 mapping and work with parallel
# ============================================================

# get command line options
# ============================================================

args=`getopt -o "r:,w:,l:" -l "resource:,wdir:,line:" -- "$@"`
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

    -l|--line)
		  line=$2
			shift 2;;

    --)
      shift
      break;;

  esac
done

# check for file availability
# ============================================================

if [[ -e $resource/resource.env && -e $resource/sbatch/hisat2.build.sbatch ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: all slurm scripts found."
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: cannot find some slurm scripts."
  exit 1
fi

if [[ -e $wdir/$line/lift.gtf && -e $wdir/$line/lift.fa ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info:starting hisat2 build"
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):info:cannot find lift gtf and fa"
  exit 1
fi

# ================
# = hisat2.build =
# ================

jobID=$(sbatch --export=env=$resource/resource.env,dir=$wdir/$line --output="$wdir"/"$line"/hisat2.build.out --error="$wdir"/"$line"/hisat2.build.err $resource/sbatch/hisat2.build.sbatch | cut -d " " -f 4)

sleep 1m
jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
while [ "$jobState" == "PD" ] || [ "$jobState" == "R" ] || [ "$jobState" == "CG" ]
do
  sleep 5m
  jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
done

reRunCounter=0
until [ `grep -i "error" "$wdir"/"$line"/hisat2.build.err | wc -l | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" "$wdir"/"$line"/hisat2.build.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 2 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting hisat2.build job."
	#sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqGRES,AllocGRES" -j $jobID >> $wdir/$line/hisat2.build.fail.out
	buildTime=$(($reRunCounter * 4 + 8))
	buildMem=$(($reRunCounter * 50 + 200))G

    rm "$wdir"/"$line"/lift.genome*

  jobID=$(sbatch --export=env=$resource/resource.env,dir=$wdir/$line --time=$buildTime:00:00 --mem=$buildMem --output="$wdir"/"$line"/hisat2.build.out --error="$wdir"/"$line"/hisat2.build.err $resource/sbatch/hisat2.build.sbatch | cut -d " " -f 4)  
	sleep 5m
	jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
	while [ "$jobState" == "PD" ] || [ "$jobState" == "R" ] || [ "$jobState" == "CG" ]
	do
		sleep 5m
		jobState=$(squeue -u tansuxu | grep $jobID | awk '{print $5}' | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `grep -i "error" "$wdir"/"$line"/hisat2.build.err | wc -l | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" "$wdir"/"$line"/hisat2.build.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed hisat2.build."
	#sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqGRES,AllocGRES" -j $jobID >> "$wdir"/"$line"/hisat2.build.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):info: error for hisat2.build."
	exit 1
fi
