
date=$(date '+%d-%b-%Y-%H-%M-%S')
mkdir -p output/$date
python simulation.py $date > output/$date/simulation_terminal_output.txt
python rates.py $date > output/$date/rates_terminal_output.txt
