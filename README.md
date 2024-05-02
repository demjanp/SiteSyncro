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
- `-max_queue_size`: Maximum queue size for parallel processing (-1 = automatic; default is -1)

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