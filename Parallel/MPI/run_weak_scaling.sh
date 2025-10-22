#!/bin/bash
#
# This script runs a WEAK SCALING test for the MPI Euler solver.
# It modifies 'N_total' in the input file based on the processor count.
#

# --- Configuration ---
EXE=./euler_mpi
PROCS=(2 4 8 16 32 64 128)
LOG_DIR="logs"
INPUT_FILE="inputs.txt"
INPUT_BAK="inputs.txt.bak"
RESULTS_FILE="weak_scaling_results.dat"

# Set the work-per-core. (e.g., 10000 cells / 16 procs = 625)
BASE_WORK_PER_CORE=625

# --- 1. Compile Code ---
echo "--- Compiling $EXE... ---"
make clean
make $EXE
if [ ! -f "$EXE" ]; then
    echo "Compilation failed. Exiting."
    exit 1
fi

# --- 2. Setup Directories and Files ---
mkdir -p $LOG_DIR
rm -f $LOG_DIR/*.log
rm -f $RESULTS_FILE
echo "--- Logs will be stored in $LOG_DIR/ ---"
echo "--- Results will be in $RESULTS_FILE ---"

# Backup the original inputs.txt
cp $INPUT_FILE $INPUT_BAK

# --- 3. Run Benchmarks ---
for P in "${PROCS[@]}"; do
    LOG_FILE="$LOG_DIR/run_np_${P}.log"
    
    # --- WEAK SCALING LOGIC ---
    # Calculate new N_total
    N_total=$(($P * $BASE_WORK_PER_CORE))
    echo "--- Modifying $INPUT_FILE: N_total = $N_total ($BASE_WORK_PER_CORE * $P) ---"
    
    # Use sed to modify the input file in-place. 
    # This finds the line starting with "N_total" and replaces the rest.
    # We use ".temp" as the backup extension to avoid conflict with our main ".bak"
    sed -i.temp "s/^\(N_total\s*\).*/\1$N_total/" $INPUT_FILE
    
    echo "--- Running with $P processes (log: $LOG_FILE)... ---"
    
    # Run the command and redirect stdout to the log file
    mpirun -np $P $EXE $INPUT_FILE > $LOG_FILE 2>&1
    
    # Clean up sed's temporary backup
    rm -f ${INPUT_FILE}.temp
done

# --- 4. Process Results ---
echo "--- Processing results... ---"

# Create the header for our data file
# Efficiency(T_2) = T_2 / T_P
echo "# Procs  N_total  TotalTime(s)  Efficiency(T_2)" > $RESULTS_FILE

T_base=0.0

for P in "${PROCS[@]}"; do
    LOG_FILE="$LOG_DIR/run_np_${P}.log"
    N_total=$(($P * $BASE_WORK_PER_CORE))
    
    # Use grep to find the "Total time" line and awk to extract the 3rd field
    T_total=$(grep "Total time:" $LOG_FILE | awk '{print $3}')
    
    if [ -z "$T_total" ]; then
        echo "Warning: Could not find 'Total time:' in $LOG_FILE"
        T_total=0
        Efficiency=0
    else
        # Get the base time (T_2) from the first run
        if [ "$P" -eq "2" ]; then
            T_base=$T_total
        fi
        
        # Use awk for robust floating-point division
        # Efficiency = T_base / T_p
        Efficiency=$(awk "BEGIN {if ($T_base > 0 && $T_total > 0) print $T_base / $T_total; else print 0;}")
    fi
    
    # Append the processed data to our results file
    echo "$P $N_total $T_total $Efficiency" >> $RESULTS_FILE
done

# --- 5. Clean up ---
echo "--- Restoring original $INPUT_FILE... ---"
mv $INPUT_BAK $INPUT_FILE

echo "Done."
echo "=============================================="
echo " Weak Scaling Results (from $RESULTS_FILE):"
echo "=============================================="
cat $RESULTS_FILE
echo "=============================================="
echo "To plot, first install Python dependencies:"
echo "  pip install matplotlib pandas"
echo "Then run:"
echo "  python3 plot_scaling.py"
