Rscript /scratch/seismo/bellinger/asteroseismology/scripts/animate.R 
#ffmpeg -start_number 260 -framerate 30 -i animate/%06d.png -c:v libx264 -r 30 animate.avi
start_number=$(ls animate | head -1 | cut -f 1 -d '.')
ffmpeg -start_number $start_number -framerate 10 -i animate/%06d.png -ab 128k -r 30 -vcodec libx264 -crf 18 -preset veryslow animate.avi


# ffmpeg -framerate 10 -i tall/%03d.png -ab 128k -r 30 -vcodec libx264 -crf 18 -preset veryslow animate.avi 



ffmpeg -start_number 003 -framerate 5 \
-i tall/%03d.png animate.avi \
-vf "tblend=average,framestep=5,setpts=10*PTS"



-filter_complex \
"[1:v][0:v]blend=all_expr='A*(if(gte(T,3),1,T/3))+B*(1-(if(gte(T,3),1,T/3)))'[v0]; \
 [2:v][1:v]blend=all_expr='A*(if(gte(T,3),1,T/3))+B*(1-(if(gte(T,3),1,T/3)))'[v1]; \
 [3:v][2:v]blend=all_expr='A*(if(gte(T,3),1,T/3))+B*(1-(if(gte(T,3),1,T/3)))'[v2]; \
 [4:v][3:v]blend=all_expr='A*(if(gte(T,3),1,T/3))+B*(1-(if(gte(T,3),1,T/3)))'[v3]; \
 [v0][v1][v2][v3]concat=n=4:v=1:a=0[v]"
