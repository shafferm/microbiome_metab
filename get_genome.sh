# Take ko table and make array of kos
string_kos="$(head -1 ko_13_5_precalculated.tab | awk '{$1=$NF=" "}1' | sed -e 's/^ *//' -e 's/ *$//')"
IFS=' ' read -a kos <<< "$string_kos"

# Take input and get line
string_counts="$(grep $1 ko_13_5_precalculated.tab | awk '{$1=$NF=" "}1' | sed -e 's/^ *//' -e 's/ *$//')"
IFS=' ' read -a counts <<< "$string_counts"

# Loop through counts and if > 0 add KO to array
genome=()
for i in "${!kos[@]}"
do
    trunc_var="$(printf %.0f ${counts[i]})"
    if [ "$trunc_var" -gt 0 ]
    then
        genome+=(${kos[i]})
    fi
done

( IFS=$'\t'; echo "${genome[*]}" )