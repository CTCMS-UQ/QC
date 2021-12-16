`QC` - A program for calculating quantum capacitance from `VASP` output.

## Usage

- The most basic invocation of the program is as follows:

```
python3 QC.py <DOS_file>
```

where you'd replace `<DOS_file>` with the filename and path of the file containing the density of states data. For example, if the DOS data 
is in the file `DOSCAR` in the same directory as the python script, you'd do:
```
    python3 QC.py DOSCAR
```
This will print the quantum capacitance (QC) at `T = 300K` for 100 evenly spaced points in the range `[-5eV,5eV]`. 

- The program has some optional arguments if you want to change the parameters of the calculation (temperature, number of points and range of potentials). Run the program with the "-h" flag to print a list of available options and how to use them:

```
python3 QC.py -h
```
For example, if you wanted to calculate the capacitance at `T=350K` you would do:

```
python3 Qc.py -T 350 
```
Everything other than the location of the data file is optional, however, and the program will run perfectly fine without them.
