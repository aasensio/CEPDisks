# CEPDisks
Set of modules to deal with disks with CEP. The routines allow the user to use a flat
disk or use disk models from external files. The disk model is placed in the
directory `disk_models`

## How to proceed

1. Generate models for each radial position running the file `generate_variable.pro` or the `generate_flat.pro`:
```
  generate_variable, disk_case, xmol
```
where `xmol` is the abundance. You probably have to edit the file to change the meaning of `disk_case`. In our
case, we set the molecular number density to $10^{-8}$ if some conditions are met.
The examples I send are for a 2-level atom (an atom with a single transition). You can change
this to a larger model atom, but it would require some modifications in the codes.

2. Solve the radiative transfer problem in each radial position by calling `calculate_variable`. 

3. The previous two steps are done for many combinations of the physical conditions using `batch_variable.pro`.
You can use the routine `batch_variable_calculate1`, for instance, as an example of how to call the previous
two routines. Additionally, you can find the routine `batch_synth1`, that solves the ray-tracing and generates
all the output to do the plots. The file `synth.pro` does this ray-tracing and some quantities are
hardwired (for example the mass of the species and the Einstein Aul coefficient). If you change the molecule
or anything else, you might have to modify them accordingly here.

4. Inside the `output_variable` directory, you find the routine `read_all.pro` that reads all the calculations
and generates and IDL save file with the line profiles.