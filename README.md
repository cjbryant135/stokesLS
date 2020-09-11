# stokesLS
A two dimesional tool for simulating the dynamics of a fluid/fluid interface in low Reynolds flow. 

## Quickstart
After downloading the source files, use 
```bash
make shearDrop_exe
```
to build all source files and an example problem. The example problem can be executed by
```bash
./shearDrop_exe <Nx> <Ny>
``` 
where *Nx*, *Ny* are the number of grid points to use in the *x* and *y* directions, respectively. 

## Authors/Acknowledgements
Colton Bryant

Level set library provided by David Chopp.
