./clean.sh
./compile.sh
./run.exe inp.2d &>results.txt
python plot_residuals.py
