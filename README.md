# PyNVEL

A Python wrapper and extension for the **USFS National Volume Estimator Library (NVEL)**.

PyNVEL provides a high-level, Pythonic interface to the low-level NVEL Fortran library via a single primary class, `VolumeCalculator`. It uses **Cython** for efficient bulk processing of tree data and for binding to the underlying NVEL routines.

In addition to the Python API, PyNVEL includes a **command-line interface (CLI)** for:

* Quick single-tree volume reports
* Exploring equation options
* Testing configuration and regional settings

> **Disclaimer**
> PyNVEL is an independent project. It is **not affiliated with or supported by** the USFS Forest Products Measurements group.

**Related resources**

* [National Volume Estimator Library overview](https://www.fs.usda.gov/managing-land/forest-management/products/measurement/volume-estimation)
* [NVEL source code](https://github.com/FMSC-Measurements/VolumeLibrary)

---

## Project Status

[![Build Python Packages](https://github.com/forest-modeling/pynvel/actions/workflows/main.yml/badge.svg)](https://github.com/forest-modeling/pynvel/actions/workflows/main.yml)

* Binary wheels available via GitHub Actions artifacts
* PyPI and Conda-Forge distribution planned

---

## Getting Started

### Installation

#### Prebuilt Wheels (Recommended)

Precompiled wheels are available as zipped artifacts from successful GitHub Actions runs.

1. Navigate to the workflow artifacts for a successful build
2. Download the archive matching your OS and Python version
3. Extract the wheel file
4. Install with `pip`

```bash
# Activate your target environment
pixi shell

# Install the extracted wheel
pip install /path/to/pynvel-*.whl
```

#### Package Managers

* **PyPI**: planned
* **Conda-Forge**: planned

---

## Usage

### API

Calculate volume for a tree.

```Python
eq = 'NVBM240202'
dbh = 18
ht = 120
vc = pynvel.VolumeCalculator(volume_eq=eq)
vc.calc(dbh_ob=dbh, total_ht=ht)

if r != 0:
    print('error code: {}'.format(r))
    print(row)

# The volume property is a dict total volume attributes
print(vc.volume['bdft_gross_prim'])
print(vc.volume['cuft_total'])

# The logs property is a list of log objects containing attributes for each log.
print(len(vc.logs))
for i,log in enumerate(vd.logs):
    print(
        f"Log: {i+1}; Diam: {round(log.small_dib, 1)}; Length {log.length};"
        f"CuFt: {round(log.cuft_gross)}; BdFt: {round(log.bdft_gross, 1)}"
        )
```

### Command-Line Interface

Display general help:

```bash
$ pynvel --help
```

#### Configuration

The configuration is stored as a JSON-formatted text file.

Export the active configuration:

```bash
$ pynvel configuration > config.json
```

Default equations and regions can be modified by editing this file.

#### Tree Volume Estimation

Generate a mini report on tree volume components.

```bash
$ pynvel volume --help
```

Estimate volume using the default equation. Species, DBH, and total height are required

```bash
$ pynvel volume -s DF -d 18 -t 120
```

```
Default species equation will be used - DF: NVBM240202
Volume Report (Version: 20231106)
---------------------------------------
Species: DF(202)
Equation: NVBM240202
DBH:             18.0
Form Class       80.0
Form Ht:          0.0
Total Ht:       120.0
Merch Ht:       101.2
Max Log:         40.0
CuFt Tot:        74.7
CuFt Merch:      70.4
BdFt Merch:     302.0
CuFt Top:         1.1
CuFt Stump:       1.8
CuFt Tip:         0.0

Product Summary
---------------
Prod Logs CuFt   BdFt    Len   Diam

Log Detail
----------
Log Prod Bole    Len     L DOB   L DIB   S DOB   S DIB   Scale   CuFt    BdFt    Int 1/4
1   1    42.0    40.0    16.3    16.3    0.0     12.0    12.0    43.6    196.0   315.0
2   1    83.0    40.0    12.0    12.0    0.0     7.6     8.0     22.7    88.0    145.0
3   2    101.0   17.0    7.6     7.6     0.0     5.0     5.0     4.1     18.0    10.0
4   3    116.0   14.0    5.0     5.0     0.0     2.0     2.0     1.1     2.0     0.0
```

Change the maximum log length to evaluate effects on merchantable volume:

```bash
$ pynvel volume -s DF -d 18 -t 120 -l 16
```

```
Default species equation will be used - DF: NVBM240202
Volume Report (Version: 20231106)
---------------------------------------
Species: DF(202)
Equation: NVBM240202
DBH:             18.0
Form Class       80.0
Form Ht:          0.0
Total Ht:       120.0
Merch Ht:       101.2
Max Log:         16.0
CuFt Tot:        74.7
CuFt Merch:      65.6
BdFt Merch:     358.0
CuFt Top:         1.1
CuFt Stump:       1.8
CuFt Tip:         0.0

Product Summary
---------------
Prod Logs CuFt   BdFt    Len   Diam

Log Detail
----------
Log Prod Bole    Len     L DOB   L DIB   S DOB   S DIB   Scale   CuFt    BdFt    Int 1/4
1   1    18.0    16.0    16.3    16.3    0.0     14.5    14.0    19.7    114.0   135.0
2   1    35.0    16.0    14.5    14.5    0.0     12.7    13.0    15.9    97.0    115.0
3   1    52.0    16.0    12.7    12.7    0.0     11.0    11.0    12.7    67.0    80.0
4   1    69.0    16.0    11.0    11.0    0.0     9.2     9.0     8.8     39.0    50.0
5   2    86.0    16.0    9.2     9.2     0.0     7.2     7.0     5.7     26.0    30.0
6   2    101.0   14.0    7.2     7.2     0.0     5.0     5.0     2.8     15.0    10.0
7   3    116.0   14.0    5.0     5.0     0.0     2.0     2.0     1.1     2.0     0.0
```

Specify an explicit NVEL equation identifier and form class as and upper stem profile point for supported profile equations:

```bash
$ pynvel volume -e F02FW3W202 -d 18 -t 120 -f 87
```

#### Stem Calculations

Estimate height to a specified **diameter inside bark (DIB)** along the bole
(Decimal values must be quoted):

```bash
$ pynvel stem-ht --help
$ pynvel stem-ht -s DF -d 18 -t 120 --stem-dib "10.4"
```

```
Stem height to 10.4" = 67.88 (54.9%)
```

Estimate **diameter inside bark** at a specified height:

```bash
$ pynvel stem-dib --help
$ pynvel stem-dib -s DF -d 18 -t 120 --stem-ht "55.5"

```

```
Stem DIB at 55.5' = 11.81
```

Both `stem-ht` and `stem-dib` accept NVEL equation arguments:

```bash
$ pynvel stem-ht -e F02FW3W202 -f 87 -d 18 -t 120 --stem-dib "10.4"
```

```
Species  is not known.
Stem height to 10.4" = 66.22 (53.4%)
```

```bash
$ pynvel stem-dib -e F02FW3W202 -f 87 -d 18 -t 120 --stem-ht "55.5"
```

```
Species  is not known.
Stem DIB at 55.5' = 11.57
```

---

## Development

### Clone the Repository

The NVEL source code is included as a Git submodule.

```bash
$ git clone --recurse-submodules https://github.com/forest-modeling/pynvel
```

For older versions of Git:

```bash
$ git clone https://github.com/forest-modeling/pynvel
$ cd pynvel
$ git submodule update --init --recursive
```

---

### Python Environment Setup

An isolated environment is strongly recommended. PyNVEL development uses **Pixi** as the preferred environment manager.

The repository includes both `pixi.toml` and `pixi.lock`.

```bash
$ cd /project/root
$ pixi install   # install dependencies
$ pixi shell     # activate environment
```

---

### Building

PyNVEL is built using:

* [Meson](https://mesonbuild.com/)
* [meson-python](https://mesonbuild.com/meson-python/)

The Cython extension and NVEL Fortran sources are compiled using GNU compilers on both Linux and Windows. All required dependencies are installed automatically when using the Pixi environment.

#### Editable Installation

An editable install automatically recompiles the extension when source files change.

> **Note**
> Recompilation will fail if the extension is currently imported in an active Python interpreter or notebook kernel.

```bash
$ python -m pip install -e . --verbose --no-build-isolation
```

#### Local Installation

```bash
$ python -m pip install . --verbose --no-build-isolation
```

#### Building Wheels

Wheels can be built for the active environment using the Python packaging build frontend:

```bash
$ python -m build --no-isolation
```

---

## Testing

A suite of unit tests is included to verify core functionality.

Tests are located in the `tests/` directory at the repository root and can be run using **pytest**:

```bash
$ cd /project/root
$ python -m pytest
```

---

## Container File

PyNVEL can be built and run as a container image. The included pynvel.containterfile provides access to 
the command line features

### Build the Image

```bash
$ podman build -t pynvel -f pynvel.containerfile
```

### Run command line tools
```bash
$ podman run --rm pynvel volume -sDF -d18 -t120
```

---

## License

See the repository license file for details.
