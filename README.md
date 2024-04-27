# SiteSyncro
Site-specific chronological modeling and synchronization

Created on 10.4.2024

<details>
<summary>Table of Contents</summary>

1. [About SiteSyncro](#about)
2. [Installation](#installation)
3. [Usage](#usage)
   * [Input File Format](#input_file)
   * [Model Class](#model_class)
   * [Sample Class](#sample_class)
4. [Developer Notes](#developer)
	* [Preparing the Virtual Environment](#venv)
	* [Building a Windows Executable](#build)
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

To run SiteSyncro on other platforms, please refer to the [Developer Notes](#developer) section on how to clone SiteSyncro from GitHub and create a virtual environment. See [Usage](#usage) section on how to run the script or import SiteSyncro as a python library.

To use SiteSyncro in your Python applications, install the latest version from the Python package index, use the following command:
```
pip install SiteSyncro
```
See [Model Class](#model_class) on usage tips.

SiteSyncro requires [OxCal](https://c14.arch.ox.ac.uk/oxcal.html) to be installed in its default location. The program is not included in the SiteSyncro package and must be installed separately. OxCal should be downloaded and installed automatically when running SiteSyncro for the first time. You can also download OxCal manually from the [OxCal website](https://c14.arch.ox.ac.uk/OxCalDistribution.zip) and unzip it in the SiteSyncro folder.

## Usage <a name="usage"></a>

### Running the Script
To use SiteSyncro, you need to run the [process.py](bin/process.py) or `sitesyncro.exe` script. This script accepts several command-line arguments. Here's a basic example of how to run the script:
```
python process.py -input data_sample.csv
```
or
```
sitesyncro.exe -input data_sample.csv
```
`process.py` & `sitesyncro.exe` accepts the following command-line arguments:
- `-h`, `--help`: Show help message and exit.
- `-directory`: Working directory for model data (default is "model").
- `-input`: The path to the input file in semicolon-separated CSV format.
- `-curve_name`: File name of the radiocarbon age calibration curve (default is "intcal20.14c").
- `-phase_model`: OxCal phase model type (can be 'sequence', 'contiguous', 'overlapping', or 'none'; default is "sequence").
- `-cluster_n`: Number of clusters to form (-1 = automatic; default is -1).
- `-by_clusters`: Flag indicating whether to update the phasing by clustering sample dates (default is 0).
- `-uniform`: Flag indicating whether to use a uniform distribution for the calendar ages (default is 0).
- `-p_value`: P-value for the randomization tests (default is 0.05).
- `-uncertainty_base`: Base uncertainty for the radiocarbon dates for the randomization tests (default is 15).
- `-npass`: Minimum number of passes for the randomization tests (default is 100).
- `-convergence`: Convergence threshold for the randomization tests (default is 0.99).
- `-max_cpus`: Maximum number of CPUs to use for parallel processing (-1 = all available; default is -1).
- `-max_queue_size`: Maximum queue size for parallel processing (default is 10000)

For example, if you want to run the script with a specific calibration curve, using uniform distributions for randomization testing and a P-value threshold of 0.01, you can do so like this:

```bash
python process.py -input data_sample.csv -curve intcal20.14c -uniform 1 -p_value 0.01
```

This will run the script with input data from data_sample.csv and use the IntCal20 calibration curve, a uniform distribution for the calendar ages, and a P-value of 0.01 for the randomization tests.

### Input File Format <a name="input_file"></a>
The input file name must be a semicolon-separated CSV file with the following 10 columns: 
1. Sample: Sample ID (required, unique identifier)
2. Context: Context ID (required)
3. Excavation Area: Excavation area ID (required)
4. C-14 Age: Radiocarbon age in years BP (required)
5. Uncertainty: Uncertainty of the radiocarbon age in C-14 years (required)
6. Phase: Phase of the sample (lower = older) (optional)
7. Earlier-Than: List of contexts that the sample is earlier (older) than (optional)
8. Long-Lived: Flag indicating whether the sample is long-lived (e.g. old wood)(required, 1 or 0)
9. Redeposited: Flag indicating whether the sample could be redeposited from a different context (required, 1 or 0)
10. Outlier: Flag indicating whether the sample is an outlier and should not be used for modeling (required, 1 or 0)

See [data_sample.csv](data_sample.csv) for an example of the input file format.

### Model Class <a name="model_class"></a>
All functions regarding modeling are encapsulated in the `Model` class. Here is a basic example of how to use it:

```python
from sitesyncro import Model

if __name__ == '__main__':
	
    # Initialize the Model object
    model = Model()
    
    # Load the data
    model.import_csv('data_sample.csv')
    
    # Process the model
    model.process()
    
    # Plot the randomization test result
    model.plot_randomized()
    
    # Plot the clustering result
    model.plot_clusters()
    
    # Save the results to a CSV file
    model.save_csv()
```
This will create the default directory `model` and generate the following files:
`model.json.gz`
`model.oxcal`
`model.js`, `model.log`, `model.txt`
`randomized.pdf`
`silhouette.pdf`
`results.csv`

#### Parameters <a name="model_parameters"></a>

The `Model` class constructor accepts the following parameters:
- `directory`: Working directory for model data (default is "model").
- `samples`: List of samples as instances of the class [Sample](#sample_class)
- `curve_name`: The name of the calibration curve to use (default is "intcal20.14c").
- `phase_model`: OxCal phase model type. Can be 'sequence', 'contiguous', 'overlapping', or 'none' (default is "sequence").
- `cluster_n`: Number of clusters to form (-1 = automatic; default is -1).
- `uniform`: Flag indicating whether to use uniform randomization (default is False).
- `p_value`: The P-value for statistical tests (default is 0.05).
- `uncertainty_base`: The base uncertainty for randomization (default is 15).
- `npass`: Minimum number of passes for the randomization tests (default is 100).
- `convergence`: Convergence threshold for the randomization tests (default is 0.99).
- `oxcal_url`: Url to download the OxCal program (default is "https://c14.arch.ox.ac.uk/OxCalDistribution.zip").

#### Methods

The `Model` class provides the following methods:
- `add_sample(sample)` or `add_sample(name, age, uncertainty, **kwargs)`: Add a sample to the model.
	- For arguments see [Sample parameters](#sample_parameters).
- `del_sample(name)`: Delete a sample from the model.
- `reset_model()`: Reset the model to the initial state.
- `save(zipped = False)`: Save the model to a file.
	- `zipped`; if True, save the model as a zipped JSON file
- `copy(directory)`: Copy the model to a new directory.
	- `directory`: Directory for the new model
- `import_csv(fname)`: Import sample data from a CSV file. See [Input File Format](#input_file) for details.
- `plot_randomized(show = False)`: Plot the randomization test results. If `show` is True, the plot is shown, otherwise it is saved as `randomized.pdf`.
- `plot_clusters(show = False)`: Plot the clustering results. If `show` is True, the plot is shown, otherwise it is saved as `silhouette.pdf`.
- `save_csv(fcsv = None)`: Save the results to a CSV file.
	- `fcsv`: File path for the CSV file. If None, `results.csv` is saved in the model directory.
- `save_outliers(fname = None)`: Saves a list of outliers to a text file
	- `fname`: File path for the text file. If None, `outliers.txt` is saved in the model directory.
- `to_oxcal()`: Save the phasing model in OxCal format as `model.oxcal`.
- `load_oxcal_data()`: Load results of OxCal modeling from `model.js`.
- `update_params(**kwargs)`: Update model parameters.
	- For keyword arguments see [Parameters](#model_parameters)
- `process_phasing(by_clusters = False)`: Update groups and phases of samples based on stratigraphic relations.
	- `by_clusters`: if True, update the phasing by clustering sample dates
- `process_outliers()`: Find dating outliers among samples which need to be removed for the model to be valid
- `process_dates()`: Calculate posteriors of sample dates based on phasing using bayesian modeling in OxCal.
- `process_randomization(max_cpus = -1, max_queue_size = 10000)`: Test if sample dates represent a uniform / normal (depending on Model.uniform parameter) distribution in time.
- `process_clustering(max_cpus = -1, max_queue_size = 10000)`: Cluster dates and using randomization testing find optimal clustering solution
- `process(by_clusters = False, max_cpus = -1, max_queue_size = 10000)`: Process the complete model
	- `by_clusters`: if True, update the phasing by clustering sample dates
	- `max_cpus`: Maximum number of CPUs to use for parallel processing (-1 = all available)
	- `max_queue_size`: Maximum queue size for parallel processing

#### Attributes

The `Model` class has the following attributes:
- `directory`: Working directory for model data
- `samples`: Dictionary of samples in format `{name: Sample, ...}`
- `curve_name`: The name of the calibration curve
- `phase_model`: OxCal phase model type
- `cluster_n`: Number of clusters to form
- `uniform`: Flag indicating whether to use uniform randomization
- `p_value`: The P-value threshold for statistical tests
- `uncertainty_base`: The base uncertainty for randomization
- `npass`: Minimum number of passes for the randomization tests
- `convergence`: Convergence threshold for the randomization tests
- `oxcal_url`: Url to download the OxCal program
- `years`: Calendar years BP corresponding to the probability distributions in format `np.array([year, ...])`
- `curve`: Calibration curve in format `np.array([[calendar year BP, C-14 year, uncertainty], ...])`, sorted by calendar years
- `uncertainties`: List of uncertainties from C-14 dates of samples
- `oxcal_data`: Results of OxCal modeling from `model.js` in format `{key: data, ...}`
- `outliers`: Returns dating outliers among samples which need to be removed for the model to be valid
- `outlier_candidates`: Returns a candidates for outliers, from which the final outliers to be eliminated were picked
	- These samples have conflicts between dating ranges and stratigraphic relationships with other samples
- `summed`: Summed probability distribution of the dating of all samples in format `np.array([p, ...])`, where p is the probability of the calendar year
- `random_p`: Calculated p-value for the randomization test
- `random_lower`: Lower bound of the randomization test in format `np.array([p, ...])`, where p is the probability of the calendar year
- `random_upper`: Upper bound of the randomization test in format `np.array([p, ...])`, where p is the probability of the calendar year
- `areas`: List of excavation areas that the samples belong to
- `contexts`: List of contexts that the samples belong to
- `groups`: Groups that the samples belong to based on stratigraphic interconnection with other samples in format `{group_name: [sample name, ...], ...}`
- `clusters`: Clusters of samples based on the similarity of their probability distributions in format `{clusters_n: [cluster: [sample name, ...], ...}, ...}`
- `cluster_means`: Mean date of the samples in each cluster in calendar years BP in format `{clusters_n: {cluster: year, ...}, ...}`
- `cluster_sils`: Silhouette score of each clustering solution in format `{clusters_n: silhouette, ...}`
- `cluster_ps`: P-value of the clustering solutions in format `{clusters_n: p, ...}`
- `cluster_opt_n`: Optimal number of clusters
- `has_data`: Boolean indicating if the model has sample data
- `is_modeled`: Boolean indicating if bayesian modeling of sample dates (process_dates) has been performed
- `is_randomized`: Boolean indicating if randomization testing of the distribution of samples (process_randomization) has been performed
- `is_clustered`: Boolean indicating if the model has been clustered (process_clustering)

### Sample Class <a name="sample_class"></a>
The `Sample` class represents a single radiocarbon sample. Here is a basic example of how to use it:

```python
from sitesyncro import Sample

if __name__ == '__main__':
	
	# Initialize the Sample object
	sample = Sample('Sample1', 1000, 50, phase = 1, earlier_than = ['Sample2'], long_lived = 1)
	
	# Print the sample data
	print(sample)
```

#### Parameters <a name="sample_parameters"></a>

The 'Sample' class constructor accepts the following parameters:
- `name`: Sample ID (required, unique identifier)
- `age`: C-14 age (years BP) for date_type 'R'; mean calendar age (years BP) for date_type 'U' (required)
- `uncertainty`: Uncertainty (years BP) for date_type 'R'; 1/2 range (years BP) for date_type 'U' (required)
- `date_type`: 'R' for radiocarbon date; 'U' for calendar date as a uniform distribution
- `long_lived`: True if sample could be older than the examined deposition event due to e.g. old wood effect
- `redeposited`: True if sample could be redeposited from a different context
- `outlier`: True if sample is an outlier and should not be used for modeling
- `context`: Name of the context where sample was found
- `area`: Excavation area
- `area_excavation_phase`: Chronological phase of the context within the excavation area (integer, lower = earlier (older) phase)
- `earlier_than`: List of names of samples which are stratigraphically later (younger) than this sample
- `curve`: Radiocarbon calibration curve in format `np.array([[calendar year BP, C-14 year, uncertainty], ...])`

#### Methods

The `Sample` class provides the following methods:
- `calibrate(curve)`: Calibrate the sample using the provided calibration curve.
	- `curve = np.array([[calendar year BP, C-14 year, uncertainty], ...])`
- `set_outlier(state)`: Set True if sample is an outlier and should not be used for modeling
- `set_group(group)`: Set the group number for the sample.
- `set_phase(phase)`: Set the phase number for the sample.
- `set_likelihood(distribution, mean = None, rng = None)`: Set the likelihood for the sample.
	- `distribution = np.array([p, ...])`
	- `rng = [from, to]`; 2-sigma (95.45%) range in calendar years BP
- `set_posterior(distribution, mean = None, rng = None, agreement = 0)`: Set the posterior for the sample.
	- `distribution = np.array([p, ...])`
	- `rng = [from, to]`; 2-sigma (95.45%) range in calendar years BP
	- `agreement`: agreement index generated by OxCal modeling
- `to_oxcal()`: Convert the sample to OxCal model format (str).
- `to_dict()`: Convert the sample data to a JSON dictionary.
- `from_dict(data)`: Load the sample data from a JSON dictionary.
- `copy()`: Create a copy of the sample instance.

#### Attributes
The `Sample` class has the following attributes:

- `name`: Sample ID (unique identifier)
- `age`: C-14 age (years BP) for date_type 'R'; mean calendar age (years BP) for date_type 'U'
- `uncertainty`: Uncertainty (years BP) for date_type 'R'; 1/2 range (years BP) for date_type 'U'
- `date_type`: 'R' for radiocarbon date; 'U' for calendar date as a uniform distribution
- `long_lived`: True if sample could be older than the examined deposition event due to e.g. old wood effect
- `redeposited`: True if sample could be redeposited from a different context
- `outlier`: True if sample is an outlier and should not be used for modeling
- `context`: Name of the context where sample was found
- `area`: Excavation area
- `area_excavation_phase`: Chronological phase of the context within the excavation area (integer, lower = earlier (older) phase)
- `earlier_than`: List of names of samples which are stratigraphically later (younger) than this sample
- `group`: Group that the sample belongs to based on stratigraphic interconnection with other samples
- `phase`: Stratigraphic phase of the sample within the group (lower = earlier (older) phase)
- `years`: Calendar years BP corresponding to the probability distributions in format `np.array([year, ...])`
- `likelihood`: Probability distribution of the dating of the sample before Bayesian modeling in format `np.array([p, ...])`, where p is the probability of the calendar year
- `posterior`: Probability distribution of the dating of the sample after Bayesian modeling in format `np.array([p, ...])`, where p is the probability of the calendar year
- `likelihood_range`: 2-sigma (95.45%) range of the dating of the sample in calendar years BP before Bayesian modeling in format `[from, to]`
- `likelihood_mean`: Mean calendar age of the sample in calendar years BP before Bayesian modeling
- `posterior_range`: 2-sigma (95.45%) range of the dating of the sample in calendar years BP after Bayesian modeling in format `[from, to]`
- `posterior_mean`: Mean calendar age of the sample in calendar years BP after Bayesian modeling
- `posterior_agreement`: Agreement index generated by OxCal modeling
- `is_calibrated`: Boolean indicating if the sample has been calibrated
- `is_modeled`: Boolean indicating if posterior has been calculated using Bayesian modeling

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

### Building a Windows Executable <a name="build"></a>
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
python -m pip install --upgrade pip
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
The executable `sitesyncro.exe` will be created in the `dist` folder.

## Contact: <a name="contact"></a>
Peter Demján (peter.demjan@gmail.com)

Institute of Archaeology of the Czech Academy of Sciences, Prague, v.v.i.

## Acknowledgements <a name="acknowledgements"></a>

Development of this software was supported by project OP JAC "Ready for the future: understanding long-term resilience of the human culture (RES-HUM)", Reg. No. CZ.02.01.01/00/22_008/0004593 of the Ministry of Education, Youth, and Sports of the Czech Republic and EU.

This software requires the [OxCal](https://c14.arch.ox.ac.uk/oxcal.html) program for bayesian modeling of the radiocarbon dates.

This software uses the following open source packages:
* [Matplotlib](https://matplotlib.org/)
* [NetworkX](https://networkx.org/)
* [NumPy](https://www.numpy.org/)
* [Requests](https://requests.readthedocs.io/)
* [Scikit-learn](https://scikit-learn.org/)
* [SciPy](https://scipy.org/)
* [tqdm](https://tqdm.github.io/)

## License <a name="license"></a>

This code is licensed under the [GNU GENERAL PUBLIC LICENSE](https://www.gnu.org/licenses/gpl-3.0.en.html) - see the [LICENSE](LICENSE) file for details