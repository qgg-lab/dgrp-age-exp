# combine mapping infomation

while read line
do
  sample=`echo $line | awk '{print $1}'`
  age=`echo $line | awk '{print $2}'`
  sex=`echo $line | awk '{print $4}'`
  
  if [[ $age == "baseline" ]] && [[ $sex == "F" ]]
  then
    cd /mnt/gs18/scratch/users/tansuxu/dgrp/lift/BaseHisat2/$sample/log
	sed 's/^ *//g' map.log|cut -d " " -f1|awk 'NR==1 {print $1}' >> /mnt/research/qgg/dgrp-age-exp/mapInfo/base.F.totalRead
	sed 's/^ *//g' map.log|cut -d " " -f1|awk 'NR==3 {print $1}' >> /mnt/research/qgg/dgrp-age-exp/mapInfo/base.F.unmap
	sed 's/^ *//g' map.log|cut -d " " -f1|awk 'NR==4 {print $1}' >> /mnt/research/qgg/dgrp-age-exp/mapInfo/base.F.unique
  elif [[ $age == "baseline" ]] && [[ $sex == "M" ]]
  then
    cd /mnt/gs18/scratch/users/tansuxu/dgrp/lift/BaseHisat2/$sample/log
	sed 's/^ *//g' map.log|cut -d " " -f1|awk 'NR==1 {print $1}' >> /mnt/research/qgg/dgrp-age-exp/mapInfo/base.M.totalRead
	sed 's/^ *//g' map.log|cut -d " " -f1|awk 'NR==3 {print $1}' >> /mnt/research/qgg/dgrp-age-exp/mapInfo/base.M.unmap
	sed 's/^ *//g' map.log|cut -d " " -f1|awk 'NR==4 {print $1}' >> /mnt/research/qgg/dgrp-age-exp/mapInfo/base.M.unique	
	
  elif [[ $age == "3wk" ]] && [[ $sex == "F" ]]
  then
    cd /mnt/gs18/scratch/users/tansuxu/dgrp/lift/3wkHisat2/$sample/log
	sed 's/^ *//g' map.log|cut -d " " -f1|awk 'NR==1 {print $1}' >> /mnt/research/qgg/dgrp-age-exp/mapInfo/3wk.F.totalRead
	sed 's/^ *//g' map.log|cut -d " " -f1|awk 'NR==3 {print $1}' >> /mnt/research/qgg/dgrp-age-exp/mapInfo/3wk.F.unmap
	sed 's/^ *//g' map.log|cut -d " " -f1|awk 'NR==4 {print $1}' >> /mnt/research/qgg/dgrp-age-exp/mapInfo/3wk.F.unique
	
  else
    cd /mnt/gs18/scratch/users/tansuxu/dgrp/lift/3wkHisat2/$sample/log
	sed 's/^ *//g' map.log|cut -d " " -f1|awk 'NR==1 {print $1}' >> /mnt/research/qgg/dgrp-age-exp/mapInfo/3wk.M.totalRead
	sed 's/^ *//g' map.log|cut -d " " -f1|awk 'NR==3 {print $1}' >> /mnt/research/qgg/dgrp-age-exp/mapInfo/3wk.M.unmap
	sed 's/^ *//g' map.log|cut -d " " -f1|awk 'NR==4 {print $1}' >> /mnt/research/qgg/dgrp-age-exp/mapInfo/3wk.M.unique
  fi

done < /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis