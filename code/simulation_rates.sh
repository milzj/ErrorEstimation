
date=$(date '+%d-%b-%Y-%H-%M-%S')
declare -a Problems=("LinearProblem" "SemilinearProblem" "BilinearProblem")

for P in "${Problems[@]}"
do
    echo $P
    mkdir -p output/$date/$P
    python simulation.py $date $P > output/$date/$P/simulation_terminal_output.txt
    python rates.py $date $P > output/$date/$P/rates_terminal_output.txt
done
