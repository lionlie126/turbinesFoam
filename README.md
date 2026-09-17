# turbinesFoam

[![DOI](https://zenodo.org/badge/4234/turbinesFoam/turbinesFoam.svg)](https://zenodo.org/badge/latestdoi/4234/turbinesFoam/turbinesFoam)
![OpenFOAM ESI v2112](https://img.shields.io/badge/OpenFOAM%20ESI-v2112-brightgreen.svg)

turbinesFoam is a library for simulating wind and marine hydrokinetic turbines
in OpenFOAM using the actuator line method.

## About this fork

This fork extends turbinesFoam with the Effective Velocity Model (EVM) for
actuator-line velocity sampling, following Schito and Zasso (2014) and Muscari
et al. (2024). It is maintained for the ESI/OpenCFD distribution of OpenFOAM;
the current development version targets OpenFOAM v2112. The fork also includes
general robustness and maintainability improvements, including preservation of
turbine azimuth across simulation restarts.

[![](https://cloud.githubusercontent.com/assets/4604869/10141523/f2e3ad9a-65da-11e5-971c-b736abd30c3b.png)](https://www.youtube.com/watch?v=THZvV4R1vow)

Be sure to check out the
[development snapshot videos on YouTube](https://www.youtube.com/playlist?list=PLOlLyh5gytG8n8D3V1lDeZ3e9fJf9ux-e).

## Installation

### Docker

Spin up an interactive shell with:

```sh
docker run --rm -it -v $PWD:/work ghcr.io/turbinesfoam/turbinesfoam
```

### Compile from source

```sh
cd $WM_PROJECT_USER_DIR
git clone https://github.com/lionlie126/turbinesFoam.git
cd turbinesFoam
./Allwmake
```

## Usage

See the tutorials located in the `tutorials` directory. The `polimi_1.2m`
tutorial provides the 1.2 m diameter PoliMi research turbine configured for
the EVM velocity-evaluation method.

### Restarting turbine simulations

Each turbine writes its cumulative azimuthal rotation to
`<time>/uniform/<turbineName>State` at normal OpenFOAM write times. When a
simulation restarts from that time directory, the turbine reads this state and
reconstructs the saved rotor position relative to its configured
`azimuthalOffset`.

For a legacy restart directory without a state file, the cumulative rotation
can be supplied manually in the turbine coefficients:

```foam
restartAngleDeg  1234.5;
```

Without saved state or `restartAngleDeg`, a nonzero-time restart emits a
warning and uses zero cumulative rotation, matching the previous behavior.

## Contributing

Pull requests are very welcome!
See the [issue tracker](https://github.com/petebachant/turbinesFoam/issues)
for more details.

## Features

`fvOptions` classes for adding actuator lines and turbines constructed from
actuator lines to any compatible solver or turbulence model, e.g.,
`simpleFoam`, `pimpleFoam`, `interFoam`, etc.

## Publications

Bachant, P., Goude, A., and Wosnik, M. (2016) [_Actuator line modeling of vertical-axis turbines_](https://arxiv.org/abs/1605.01449). arXiv preprint 1605.01449.

Schito, P., and Zasso, A. (2014) [_Actuator forces in CFD: RANS and LES modeling in OpenFOAM_](https://doi.org/10.1088/1742-6596/524/1/012160). Journal of Physics: Conference Series, 524, 012160.

Muscari, C., Schito, P., Viré, A., Zasso, A., and van Wingerden, J.-W. (2024) [_The effective velocity model: An improved approach to velocity sampling in actuator line models_](https://doi.org/10.1002/we.2894). Wind Energy, 27(5), 447-462.

## How to cite

The latest release of turbinesFoam can be cited via DOI thanks to Zenodo: [![DOI](https://zenodo.org/badge/4234/turbinesFoam/turbinesFoam.svg)](https://zenodo.org/badge/latestdoi/4234/turbinesFoam/turbinesFoam)

## Acknowledgements

This work was funded through a National Science Foundation CAREER award,
principal investigator Martin Wosnik ([NSF CBET
1150797](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1150797), Energy for
Sustainability, original program manager Geoffrey A. Prentice, current program
manager Gregory L. Rorrer).

OpenFOAM is free, open source software for computational fluid dynamics (CFD),
developed primarily by [CFD Direct](http://cfd.direct), on behalf of the
[OpenFOAM](http://openfoam.org) Foundation.

Interpolation, Gaussian projection, and vector rotation functions adapted from
NREL's [SOWFA](https://github.com/NREL/SOWFA).
