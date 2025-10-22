#!/bin/bash
#
# This script runs a STRONG SCALING test for the MPI Euler solver.
# It uses a FIXED problem size (from inputs.txt) and varies the core count.
#

# --- Configuration ---
EXE=./euler_mpi
PROCS=(2 4 8 16 32 64 128)
LOG_DIR="logs_strong"
INPUT_FILE="inputs.txt"
RESULTS_FILE="strong_scaling_results.dat"

# --- 1. Compile Code ---
echo "--- Compiling $EXE... ---"
# Make sure we have the latest executable
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
echo "--- Using fixed problem size from $INPUT_FILE ---"

# --- 3. Run Benchmarks ---
for P in "${PROCS[@]}"; do
    LOG_FILE="$LOG_DIR/run_np_${P}.log"
    
    echo "--- Running with $P processes (log: $LOG_FILE)... ---"
    
    # Run the command and redirect stdout to the log file
    mpirun -np $P $EXE $INPUT_FILE > $LOG_FILE 2>&1
done

# --- 4. Process Results ---
echo "--- Processing results... ---"

# Create the header for our data file
# Speedup(T_2) = T_2 / T_P
echo "# Procs  TotalTime(s)  Speedup(T_2)" > $RESULTS_FILE

T_base=0.0
BASE_PROCS=${PROCS[0]} # Get the first element (e.g., 2)

for P in "${PROCS[@]}"; do
    LOG_FILE="$LOG_DIR/run_np_${P}.log"
    
    # Use grep to find the "Total time" line and awk to extract the 3rd field
    T_total=$(grep "Total time:" $LOG_FILE | awk '{print $3}')
    
    if [ -z "$T_total" ]; then
        echo "Warning: Could not find 'Total time:' in $LOG_FILE"
        T_total=0
        Speedup=0
    else
        # Get the base time (T_base) from the first run
        if [ "$P" -eq "$BASE_PROCS" ]; then
            T_base=$T_total
        fi
        
        # Use awk for robust floating-point division
        # Speedup = T_base / T_p
        Speedup=$(awk "BEGIN {if ($T_base > 0 && $T_total > 0) print $T_base / $T_total; else print 0;}")
    fi
    
    # Append the processed data to our results file
    echo "$P $T_total $Speedup" >> $RESULTS_FILE
done

# --- 5. Clean up ---
echo "Done."
echo "=============================================="
echo " Strong Scaling Results (from $RESULTS_FILE):"
echo "=============================================="
cat $RESULTS_FILE
echo "=============================================="
echo "To plot, first install Python dependencies:"
echo "  pip install matplotlib pandas"
echo "Then run:"
echo "  python3 plot_strong_scaling.py"