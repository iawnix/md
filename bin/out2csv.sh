FILE_=$1

grep -v "#" ./${FILE_}.out | awk '{printf "%s,%s\n",$1,$2}' > ${FILE_}.csv 
