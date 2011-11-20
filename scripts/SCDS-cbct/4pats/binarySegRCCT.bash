#seeds
#pat21:  61 26 17   199 126  14
#pat25:  84 124  11   191 124 13    # uppder lung tumor
#pat31:   56 126  75   193  136  54
#pat20: 88 126 12  216 126 17


## threshold

for j in  20  25 31
do
	cd /stage/sharonxx/proj/mskcc/Patient$j/rcct/image
	mkdir binary-inter
	for i in 00 10 20 30 40 50 60 70 80 90 
	do 
		ThresholdImage 3 ./gray-inter/cinePhase$i.mhd ./binary-inter/lung-bin-cinePhase$i.mhd  -900 -300


	done
done


## connected regions
for j in 20
do
cd /stage/sharonxx/proj/mskcc/Patient$j/rcct/image
for i in 00 10 20 30 40 50 60 70 80 90 
do 
ConnectedThreshold ./binary-inter/lung-bin-cinePhase$i.mhd  ./binary-inter/lung-bin-cinePhase$i.mhd   1  2 88 126 12  216 126 17
done
done





for j in   21 
do
cd /stage/sharonxx/proj/mskcc/Patient$j/rcct/image
mkdir binary-inter
for i in 00 10 20 30 40 50 60 70 80 90 
do 
ThresholdImage 3 ./gray-inter/cinePhase$i.mhd ./binary-inter/lung-bin-cinePhase$i.mhd  -1000 -300
done
done


for j in 21
do
cd /stage/sharonxx/proj/mskcc/Patient$j/rcct/image
for i in 00 10 20 30 40 50 60 70 80 90 
do 
ConnectedThreshold ./binary-inter/lung-bin-cinePhase$i.mhd  ./binary-inter/lung-bin-cinePhase$i.mhd   1  2  199 126  14
done
done




for j in 25
do
	cd /stage/sharonxx/proj/mskcc/Patient$j/rcct/image
	for i in 00 10 20 30 40 50 60 70 80 90 
	do 
		ConnectedThreshold ./binary-inter/lung-bin-cinePhase$i.mhd  ./binary-inter/lung-bin-cinePhase$i.mhd   1  2  84 124  11   191 124 13

	done
done



for j in 31
do
	cd /stage/sharonxx/proj/mskcc/Patient$j/rcct/image
	for i in 00 10 20 30 40 50 60 70 80 90 
	do 
		ConnectedThreshold ./binary-inter/lung-bin-cinePhase$i.mhd  ./binary-inter/lung-bin-cinePhase$i.mhd   1  2  56 126  75   193  136  54

	done
done





