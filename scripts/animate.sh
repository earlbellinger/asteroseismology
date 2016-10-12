Rscript /scratch/seismo/bellinger/asteroseismology/scripts/animate.R 
#ffmpeg -start_number 260 -framerate 30 -i animate/%06d.png -c:v libx264 -r 30 animate.avi
start_number=$(ls animate | head -1 | cut -f 1 -d '.')
ffmpeg -start_number $start_number -framerate 30 -i animate/%06d.png -ab 128k -r 30 -vcodec libxvid animate.avi
