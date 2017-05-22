maybe_sub.sh -p 16 ./dispatch.sh -M 1.2 -Y 0.294754 -Z 0.0176894 -d 0.01 -D 0 -d overshoot
maybe_sub.sh -p 16 ./dispatch.sh -M 1.2 -Y 0.294754 -Z 0.0176894 -o 0    -D 0 -d novershoot2

maybe_sub.sh -p 16 ./dispatch.sh -M 1.15 -Y 0.294754 -Z 0.0176894 -o 0    -D 0 -d lowerM

maybe_sub.sh -p 16 ./dispatch.sh -D 0 -d sun

maybe_sub.sh -p 16 ./dispatch.sh -M 1.2 -o 0 -D 0 -d solar1_2

maybe_sub.sh -p 16 ./dispatch.sh -M 1.11 -Y 0.293544 -Z 0.0176313 -o 0 -D 0 -d gemma2
maybe_sub.sh -p 16 ./dispatch.sh -M 1.11 -o 0 -D 0 -d gemma0

#maybe_sub.sh -p 8 ./dispatch.sh -M 1.2 -Y 0.294754 -Z 0.0176894 -D 0 -s 1 -d step 
#maybe_sub.sh -p 8 ./dispatch.sh -M 1.2 -Y 0.294754 -Z 0.0176894 -D 0 -s 0 -d exp 

