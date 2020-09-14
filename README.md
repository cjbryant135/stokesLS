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
where *Nx*, *Ny* are the number of grid points to use in the *x* and *y* directions, respectively. On the first run, the user will be prompted to create a parameter file specifying the initial position/radius of the drop. Answer `y` and input drop position/radius. Note that the default domain size is `[-1,1]x[-1,1]`.

## Details
See docs/stokesLS.pdf for mathematical details. Code documentation can be found [here](https://cjbryant135.github.io/stokesLS/html/md_README.html).

## Authors/Acknowledgements
Colton Bryant

Level set library provided by David Chopp.
