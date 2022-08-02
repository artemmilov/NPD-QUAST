# About NPD-Quast
NPD-Quast is designed for analyzing natural peptide identificator\`s reports of the provided data and comparing them to each other.

* Copability: NPD-Quast is Python 3.9+ compliant
* Prerequisites: Apart from a working version of Python, the following non-standard libraries need to be installed: conda, numpy and matplotlib. Also you need to install RDKit if you want to run NPD-Quast whithout `npd-quast-env` activated.
* Disclaimer: NPD-Quast should work on all platforms, but was only tested in a Linux (Ubuntu 20.04) system.

## How to instal NPD-Quast?
Just clone this repository on your local machine and fill `npd_quast.ini` file. Note that you need to write location just that tools, which you want to run.

## How to construct work folder?
The work folder should contain subdir with challenges, subdir with reports and file with true answers. There are three examples of work folders in `sample` dir.

### What are challenges are?
Each challenge consist of database and spectra list. The tool under test makes identifications between spectra and molecules from database, corresponding to this challenge.

### Subdir with challenges
Firstly you need to make subsir for every challenge you have. In each challenge dir write the database file in format `csv` and make subdir named `spectra`. You have to make `mgf` file for every spectrum you have for this challenge and put them to the `spectra` subdir. Note each spectrum file contain only one spectrum.

### Subdir with reports
Just make empty subdir named `reports`.

### True answers
This file should be named "true.answers.txt". Each line of it corresponds to the true identification. All true identification should be written there. So if the `spectrum_name` spectrum from `challenge_name` challenge matches to the molecule with `true_answer_inchi_key` inchi key you need to write line `challenge_name'\t'spectrum_name'\t'true_answer_smiles'\t'true_answer_inchi` in this file. Note that you need only first 11 symbols from inchi key.

## How to run NPD-Quast?
To run NPD-Quast you firstly need to construct work folder which structure is written above. Also, you need all dependencies to be installed. We recommend you activate a `npd-quast-env.yml` environment which is located in `envs` dir. So, just run this command:

```conda activate npd-quast-env```

### Running supported tool
If you want to run your data by one of the supported tool (assume it\`s named `tool_name`) run such command (here we decided to name page with report `report_name`, you can name it as you want):

```python npd_quast.py run_n_report tool_name report_name path_to_work_folder```

### Changing tool parameters
To change tool parameters you need to copy json file with original params and change it as you need. After run this command:

```python npd_quast.py run_n_report tool_name report_name path_to_work_folder --config path_to_your_configuration```

### Adding unsupported tool reports
If you want to report unsupported tool you firstly need to run it on wour machine by yourself. After that make subdir `report_name` in the `reports` dir and write there `tool_answers.txt` file. It structure is similar as `true_answers.txt` structure. 

Each line of it corresponds to the tool identification. All tool identification should be written there. So if the `spectrum_name` spectrum from `challenge_name` challenge matches to the molecule with `tool_answer_inchi_key` inchi key you need to write line `challenge_name'\t'spectrum_name'\t'tool_answer_inchi'\t'match_score` in this file. Note that you need only first 11 symbols from inchi key.

After that just run command:

```python npd_quast.py compile_reports path_to_work_folder```

## Supported tools
* Dereplicator (https://github.com/ablab/npdtools)
* Dereplicator+ (https://github.com/ablab/npdtools)
* MAGMa+ (https://github.com/savantas/MAGMa-plus)
* Sirius (https://github.com/boecker-lab/sirius)
