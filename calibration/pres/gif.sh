#rm -rf *.png
Rscript hrs.R
convert -delay 6 -loop 0 HR/*.png HR/HR.gif
convert -loop 0 HR/HR.gif -delay 200 HR/30HR.png HR_pause.gif

convert -delay 6 -loop 0 u/*.png u/u.gif
convert -loop 0 u/u.gif -delay 200 u/30u.png u_pause.gif


