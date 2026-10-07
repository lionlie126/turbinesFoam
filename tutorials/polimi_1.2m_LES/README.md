# PoliMi 1.2 m turbine — LES case

This case provides a production-scale large-eddy simulation setup for the
1.2 m diameter, three-bladed PoliMi research turbine. It complements the
smaller `polimi_1.2m` tutorial and retains the refined wake mesh and sampling
configuration used for LES analysis.

The operating point and actuator-line settings are:

- free-stream velocity: 4 m/s;
- tip-speed ratio: 7.5;
- 0.6 m rotor radius and 38 actuator elements per blade;
- EVM velocity evaluation with five sampling points;
- isotropic Gaussian force projection with a chord-based radius;
- dynamic stall, filtered lifting-line corrections and end effects disabled.

The LES setup uses the WALE model with a `cubeRootVol` filter width, a
second-order backward time scheme and LUST convection for velocity. The mesh
contains a level-2 cylindrical near-wake region and level-3 refinement around
the supplied `shearLayer.stl` surface. The aerodynamic data comprise 39
spanwise multi-Reynolds-number lift and drag tables.

The meshing stage uses the 10-domain decomposition in
`system/decomposeParDict.snappy`. The parallel flow solution uses the 40-domain
decomposition in `system/decomposeParDict`.

Run a serial flow solution after parallel meshing with:

```bash
./Allrun
```

Run both meshing and the flow solution in parallel with:

```bash
./Allrun -parallel
```

The case is considerably more expensive than the compact `polimi_1.2m`
tutorial and is intended for a multi-core workstation or cluster. Adjust both
decomposition dictionaries when using different processor counts.

After a decomposed run, the wake-line data can be sampled with:

```bash
mpirun -np 40 postProcess -parallel \
    -dict system/sampleLinesDict -latestTime
```

The case targets the OpenCFD/ESI OpenFOAM line, with v2112 as the maintained
reference version.
