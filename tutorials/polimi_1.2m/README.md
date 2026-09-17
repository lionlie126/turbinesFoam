# PoliMi 1.2 m turbine

This case provides the 1.2 m diameter, three-bladed PoliMi research turbine as
an axial-flow actuator-line example. The turbine has a 0.6 m rotor radius, 38
actuator elements per blade, and 39 spanwise aerodynamic profiles.

The supplied operating point follows the current PoliMi configuration:

- free-stream velocity: 4.05 m/s;
- tip-speed ratio: 7.5;
- EVM velocity evaluation with five points;
- isotropic Gaussian force projection with a grid-based radius;
- dynamic stall, filtered lifting-line corrections, and end effects disabled.

The lift and drag tables cover chord Reynolds numbers from 30,000 to 250,000.
The source data identify the full-angle-of-attack polars as derived from
Sheldahl and Klimas (1981). Pitching-moment tables are omitted because the
available tables contain zero values and the model treats omitted moment data
as zero.

The computational domain and runtime are intentionally smaller than the
research production case. They are intended to demonstrate configuration of
the turbine and EVM implementation, not to reproduce a validated performance
curve or a mesh-convergence study.

Run the case in the same way as the other tutorials:

```bash
./Allrun
```

For a parallel run:

```bash
./Allrun -parallel
```

The case targets the OpenCFD/ESI OpenFOAM line, with v2112 as the maintained
reference version.
