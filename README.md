# SiteSyncro
Site-specific chronological modeling and synchronization

Created on 10.4.2024

<details>
<summary>Table of Contents</summary>

1. [About SiteSyncro](#about)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Developer Notes](#developer)
5. [Contact](#contact)
6. [Acknowledgements](#acknowledgements)
7. [License](#license)

</details>

## About SiteSyncro <a name="about"></a>
SiteSyncro is a Python-based tool designed for site-specific chronological modeling and synchronization based on radiocarbon dates from stratigraphically and/or spatially linked archaeological contexts. It uses the [OxCal](https://c14.arch.ox.ac.uk/oxcal.html) program for bayesian modelling and a method of temporal clustering of the modelled C-14 dates to determine whether they represent separate events, or phases in time.
 
The input data represent radiocarbon dates and their stratigraphic relations. Here is a brief overview of the processing workflow:
1. Bayesian modeling of distributions of the C-14 dates based on stratigraphic constrains using the [OxCal](https://c14.arch.ox.ac.uk/oxcal.html) program.
2. Randomization testing of the null hypothesis that the observed C-14 dates represent a normal / uniform distribution
3. Temporal clustering of the modelled C-14 dates. All possible clusterings are produced and tested for randomness. The optimal number of clusters is selected for further analyses.
4. Updating of the bayesian model phasing based on the temporal clustering.

For a detailed overview of the clustering method see:

Demján, P., & Pavúk, P. (2021). CLUSTERING OF CALIBRATED RADIOCARBON DATES: SITE-SPECIFIC CHRONOLOGICAL SEQUENCES IDENTIFIED BY DENSE RADIOCARBON SAMPLING. Radiocarbon, 63(2), 429-438. [doi:10.1017/RDC.2020.129](https://doi.org/10.1017/RDC.2020.129)


## Installation <a name="installation"></a>

[Windows executable](https://github.com/demjanp/SiteSyncro/releases/latest) is available for users who do not want to install Python and dependencies.

To run SiteSyncro on other platforms, please refer to the [Developer Notes](#developer) section on how to clone SiteSyncro from GitHub and create a virtual environment. See [Usage](#usage) section on how to run the script.

To use SiteSyncro in your Python applications, install the latest version from the Python package index, use the following command:
```
pip install SiteSyncro
```
See [SiteSyncro Class](#sitesync_class) on usage tips.

SiteSyncro requires [OxCal](https://c14.arch.ox.ac.uk/oxcal.html) to be installed in its default location. The program is not included in the SiteSyncro package and must be installed separately. OxCal should be downloaded and installed automatically when running SiteSyncro for the first time. You can also download OxCal manually from the [OxCal website](https://c14.arch.ox.ac.uk/OxCalDistribution.zip) and unzip it in the SiteSyncro folder.

## Usage <a name="usage"></a>
To use SiteSyncro, you need to run the [process.py](bin/process.py) or `sitesyncro.exe` script. This script accepts several command-line arguments. Here's a basic example of how to run the script:
```
python process.py data_sample.csv
```
or
```
sitesyncro.exe data_sample.csv
```
`process.py` & `sitesyncro.exe` accepts the following command-line arguments:

- `input`: The path to the input file in semicolon-separated CSV format (required).
- `-result`: The directory path to store the results (default is "result").
- `-existing`: Flag indicating whether to use existing results (default is 0).
- `-curve`: File name of the radiocarbon age calibration curve (default is "intcal20.14c").
- `-model`: OxCal model type (can be 'sequence', 'contiguous', 'overlapping', or 'none'; default is "sequence").
- `-n`: Number of clusters to form (-1 = automatic; default is -1).
- `-uniform`: Flag indicating whether to use a uniform distribution for the calendar ages (default is 0).
- `-p_value`: P-value for the randomization tests (default is 0.05).
- `-uncert_base`: Base uncertainty for the radiocarbon dates for the randomization tests (default is 15).
- `-npass`: Minimum number of passes for the randomization tests (default is 100).
- `-convergence`: Convergence threshold for the randomization tests (default is 0.99).
- `-max_cpus`: Maximum number of CPUs to use for parallel processing (-1 = all available; default is -1).
- `-max_queue_size`: Maximum queue size for parallel processing (default is 100).

For example, if you want to run the script with a specific calibration curve, using uniform distributions for randomization testing and a P-value threshold of 0.01, you can do so like this:

```bash
python process.py data_sample.csv -curve intcal20.14c -uniform 1 -p_value 0.01
```

This will run the script with input data from data_sample.csv and use the IntCal20 calibration curve, a uniform distribution for the calendar ages, and a P-value of 0.01 for the randomization tests.

The input file name must be a semicolon-separated CSV file with the following 8 columns: 
1. Sample: Sample ID (required, unique identifier)
2. Context: Context ID (required)
3. Excavation Area: Excavation area ID (required)
4. C-14 Age: Radiocarbon age in years BP (required)
5. Uncertainty: Uncertainty of the radiocarbon age in C-14 years (required)
6. Phase: Phase of the sample (lower = older) (optional)
7. Earlier-Than: List of contexts that the sample is earlier (older) than (optional)
8. Long-Lived / Redeposited: Flag indicating whether the sample is long-lived (e.g. old wood) or redeposited (e.g. from an older context) (required, 1 or 0)

See [data_sample.csv](data_sample.csv) for an example of the input file format.

## Developer Notes <a name="developer"></a>
### Preparing the Virtual Environment <a name="venv"></a>
SiteSyncro requires [Python 3.10](https://www.python.org/downloads/).

To prepare a Python virtual environment open a terminal or command prompt window and type the following commands:

```
git clone https://github.com/demjanp/SiteSyncro.git

cd sitesyncro

python -m venv .venv

.venv\Scripts\activate.bat

python.exe -m pip install --upgrade pip

pip install -e .
```

See [Usage](#usage) on further instructions how to run the script.

### Building a Windows Executable
To build a Windows executable, open a terminal or command prompt window and change to the `sitesyncro` folder: 
```
cd SiteSyncro\sitesyncro
```
Activate the virtual environment:
```
.venv\Scripts\activate.bat
```
Then type the following commands (this only has to be done once per virtual environment):
```
python.exe -m pip install --upgrade pip
python -m pip install --upgrade build
pip install twine
pip install pyinstaller==6.6.0
```
To build the executable, type:
```
pip install -e .
python -m build
pyinstaller sitesyncro.spec
```
You will find the executable `sitesyncro.exe` in the `dist` folder.

### SiteSyncro Class <a name="sitesync_class"></a>
All functions of SiteSyncro are encapsulated in the `SiteSyncro` class. Here is a basic example of how to use it:

```python
from sitesyncro import SiteSyncro

if __name__ == '__main__':
	
    # Initialize the SiteSyncro object
    ssync = SiteSyncro(input='data_sample.csv')
    
    # Load the data
    ssync.load_data()
    
    # Process the model
    ssync.process()
    
    # Plot the randomization test result
    ssync.plot_randomized()
    
    # Plot the clustering result
    ssync.plot_clusters()
    
    # Save the results to a CSV file
    ssync.save_results_csv()
```

#### Parameters <a name="parameters"></a>

The `SiteSyncro` class accepts the following parameters:
- `input`: The path to the input file.
- `result`: The name of the result directory (default is "result").
- `existing`: Boolean indicating whether to use existing results (default is False).
- `curve`: The name of the calibration curve to use (default is "intcal20.14c").
- `model`: The name of the model to use. Can be 'sequence', 'contiguous', 'overlapping', or 'none' (default is "sequence").
- `n`: The number of clusters (-1 = automatic) (default is -1).
- `uniform`: Boolean indicating whether to use uniform randomization (default is False).
- `p_value`: The P-value for statistical tests (default is 0.05).
- `uncert_base`: The base uncertainty for randomization (default is 15).
- `npass`: The number of passes for the clustering algorithm (default is 100).
- `convergence`: The convergence criterion for the clustering algorithm (default is 0.99).
- `max_cpus`: The maximum number of CPUs to use (default is -1).
- `max_queue_size`: The maximum size of the queue for parallel processing (default is 100).

#### Methods

The `SiteSyncro` class provides the following methods:
- `load_data(**kwargs)`: Loads data from the input file and updates the internal data. The input file must be a CSV file with the following columns: Context, Excavation Area, C-14 Age, Uncertainty, Phase (lower = older), Earlier-Than, Long-Lived / Redeposited. For possible keyword arguments see [Parameters](#parameters).
- `process(**kwargs)`: Processes the data. For possible keyword arguments see [Parameters](#parameters).
- `plot_randomized(show=False)`: Plots the randomized data. If `show` is True, the plot is shown. If False, the plot is saved to a file.
- `plot_clusters(show=False)`: Plots the clustering data. If `show` is True, the plot is shown. If False, the plot is saved to a file.
- `save_results_csv(fcsv=None)`: Saves the results to a CSV file. If `fcsv` is None, a default file path is used.

#### Attributes

The `SiteSyncro` class has the following attributes:
- `oc_data` (OxCalData): Used to store and manage the data after bayesian modelling of the radiocarbon dates with the OxCal program.
- `rnd_data` (RandomizeData): Used to store and manage the data related to the randomization testing of the null hypothesis that the observed C-14 dates represent a normal/uniform distribution.
- `clu_data` (ClusterData): Used to store and manage the data related to the temporal clustering of the modelled C-14 dates.
- `oc_clu_data` (OxCalData): Used to store and manage the modelled data updated according to the temporal clustering of the C-14 dates.

All variables generated by the OxCal program can be accessed using `oc_data[key]`. Use `print(oc_data)` to see a list of available keys.

`OxCalData` class methods:
- `__init__(fname=None)`: If a filename of an OxCal .js file is provided, it reads and stores the data.
- `has_data(self)`: Returns True if the data dictionary is not empty, False otherwise.
- `set_priors(samples, context_samples, context_area, areas, groups, phases, earlier_than, long_lived, r_dates, dates, model, curve)`: Sets the priors for the bayesian modeling.
    ```
    samples = [sample, ...]
    context_samples = {context: [sample, ...], ...}
    context_area = {context: area, ...}
    areas = [area, ...]
    groups = {group: [sample, ...], ...}
    phases = {group: {sample: phase, ...}, ...}
    long_lived = {sample: True/False, ...}
    r_dates = {sample: (age, uncertainty), ...}
    dates = {sample: (age, uncertainty, 'R'), sample: (CE from, CE to, 'U'), ...}; R = radiocarbon, U = uniform
    ```
- `get_likelihoods()`: Returns the likelihood distributions for each sample.
    Format: `{name: {'prob': [p, ...], 'years': [calendar year CE, ...], 'mean': mean, 'range': [min year CE, max year CE]}}`
    Ranges are calculated at the 95.45% (2-sigma) level.
- `get_posteriors()`: Returns the posterior distributions for each sample. 
    Format: `{name: {'prob': [p, ...], 'years': [calendar year CE, ...], 'mean': mean, 'range': [min year CE, max year CE], 'agreement': agreement index}}`	
    Ranges are calculated at the 95.45% (2-sigma) level.
- `get_samples()`: Returns the sample names.
    Format: `[name, ...]`
- `get_context_samples()`: Returns the contexts and their contained samples.
    Format: `{context: [sample, ...], ...}`
- `get_context_area()`: Returns the contexts and their areas.
    Format: `{context: area, ...}`
- `get_areas()`: Returns the areas.
    Format: `[area, ...]`
- `get_groups()`: Returns the groups based on stratigraphic links.
    Format: `{group: [sample, ...], ...}`
- `get_phases()`: Returns the phases for each group and sample. Lower = earlier (older) phase.
    Format: `{group: {sample: phase, ...}, ...}`
- `get_earlier_than()`: Returns the earlier-than relations between samples.
    Format: `matrix[n_samples x n_samples] = [True/False, ...]`. Sample in row is earlier than sample in column based on stratigraphy.
- `get_long_lived()`: Returns the long-lived value for each sample.
    Format: `{sample: True/False, ...}`
- `get_r_dates()`: Returns the radiocarbon dates for each sample.
    Format: `{sample: (age, uncertainty), ...}`
- `get_dates()`: Returns the date for bayesian modeling for each sample. Long-lived or redeposited samples are represented by a uniform distribution where the lower boundary is the lower boundary of the 2-sigma range of the calibrated date and the upper boundary is 1950 CE. 
    Format: `{sample: (age, uncertainty, 'R'), sample: (CE from, CE to, 'U'), ...}`; R = radiocarbon, U = uniform
- `get_model()`: Returns the model type used for the bayesian modeling.
    Values: 'sequence', 'contiguous', 'overlapping', 'none'
- `get_curve()`: Returns the calibration curve used.
    Format: curve name (see OxCal\bin folder for available curves)

`ClusterData` class methods:
- `__init__(clusters = {}, means = {}, sils = {}, ps = {}, p_value = None, opt_n = None)`
    ```
    means = {n: mean, ...}
    sils = {sample: sil, ...}
    ps = {n: p, ...}
    p_value = p-value threshold for randomization testing
    opt_n = optimal number of clusters
	```	
- `has_data()`: Returns True if the clusters dictionary is not empty, False otherwise.
- `get_clusters()`: Returns the clusters for each number of clusters.
    Format: `{n: {sample: cluster, ...}, ...}`
- `get_means()`: Returns the mean silhouette score for each number of clusters.
    Format: `{n: mean, ...}`
- `get_sils()`: Returns the silhouette score for each sample.
    Format: `{sample: sil, ...}`
- `get_ps()`: Returns the p-value for each number of clusters.
    Format: `{n: p, ...}`
- `get_p_value()`: Returns the p-value threshold for randomization testing.
- `get_opt_n()`: Returns the optimal number of clusters.

`RandomizeData` class methods:
- `__init__(years = None, sum_obs = None, uniform = None, p = None, p_value = None, sums_rnd_lower = None, sums_rnd_upper = None)`
    ```
    years = [year, ...]. The years are in calendar years BP.
    sum_obs = [p, ...] in order of the years.
    uniform = True/False
    p = calculated p-value
    p_value = p-value threshold
    sums_rnd_lower = [p, ...] in order of the years.
    sums_rnd_upper = [p, ...] in order of the years.
	```
- `has_data()`: Returns True if the years attribute is not None, False otherwise.
- `get_years()`: Returns the years for the randomization test.
    Format: `[year, ...]`. The years are in calendar years BP.
- `get_sum_obs()`: Returns the sum of the observed distributions.
    Format: `[p, ...]` in order of the years.
- `get_uniform()`: Returns True if the uniform distribution was used.
- `get_p()`: Returns the calculated p-value for the randomization test.
- `get_p_value()`: Returns the p-value threshold for the randomization test.
- `get_sums_rnd_lower()`: Returns the lower boundary of the randomized distributions.
    Format: `[p, ...]` in order of the years.
- `get_sums_rnd_upper()`: Returns the upper boundary of the randomized distributions.
    Format: `[p, ...]` in order of the years.

## Contact: <a name="contact"></a>
Peter Demján (peter.demjan@gmail.com)

Institute of Archaeology of the Czech Academy of Sciences, Prague, v.v.i.

## Acknowledgements <a name="acknowledgements"></a>

Development of this software was supported by project OP JAC "Ready for the future: understanding long-term resilience of the human culture (RES-HUM)", Reg. No. CZ.02.01.01/00/22_008/0004593 of the Ministry of Education, Youth, and Sports of the Czech Republic and EU.

This software requires the [OxCal](https://c14.arch.ox.ac.uk/oxcal.html) program for bayesian modeling of the radiocarbon dates.

This software uses the following open source packages:
* [NumPy](https://www.numpy.org/)
* [SciPy](https://scipy.org/)
* [Matplotlib](https://matplotlib.org/)
* [Scikit-learn](https://scikit-learn.org/)
* [tqdm](https://tqdm.github.io/)
* [Requests](https://requests.readthedocs.io/)
* [NetworkX](https://networkx.org/)

## License <a name="license"></a>

This code is licensed under the [GNU GENERAL PUBLIC LICENSE](https://www.gnu.org/licenses/gpl-3.0.en.html) - see the [LICENSE](LICENSE) file for details