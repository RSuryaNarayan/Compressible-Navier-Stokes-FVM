make clean && make && ./gridGenRamp inp.ramp && ./readGrid inp.ramp &> results.txt && python plotGridCells.py inp.ramp && python plotSolution.py inp.ramp
