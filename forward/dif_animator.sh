gifs=$(find tall -name nuclear-*.png|sort|xargs)
_CONVERT="convert "

for gif in $gifs; do
    name=$(echo ${gif##*/})
    delay=${name:8} #6}
    delay=${delay::-4}
    #delay=$(echo "$delay / 10" | bc -l)
    if (( $(bc <<< "$delay > 0.01") )); then 
        _CONVERT="${_CONVERT} -delay $delay $gif " 
    fi
done;

_CONVERT="${_CONVERT} -loop 0 -layers Optimize diffusion-nuclear.gif"

eval $_CONVERT

