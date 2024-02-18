
date=$(date '+%d-%b-%Y-%H-%M-%S')
mkdir -p output/$date
python simulation.py $date 2>&1 | tee output/$date/simulation_terminal_output.txt
python rates.py $date 2>&1 | tee output/$date/rates_terminal_output.txt
