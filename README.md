# NPD-Quast
NPD-Quast is designed for analyzing natural peptide identificator\`s reports of the provided data and comparing them to each other.

* Copability: NPD-Quast is Python 3.9+ compliant
* Prerequisites: Apart from a working version of Python, the following non-standard libraries need to be installed: conda, numpy and matplotlib. Also you need to install RDKit if you want to run NPD-Quast whithout `npd-quast-env` activated.
* Disclaimer: NPD-Quast should work on all platforms, but was only tested in a Linux (Ubuntu 20.04) system.

## Installing NPD-Quast
Just clone this repository on your local machine.

## Running NPD-Quast
* Fill `npd_quast.ini` file. Note that you need to write location just that tools, which you want to run.
* We recommend you to activate `npd-quast-env`. It is located in "envs" dir.
* To add new tool report type such command in terminal  
```
	python <npdquast_installation_dir>/npd_quast.py new_report <tool_name> <work_folder>
```
* To add total report of existed tool reports in work folder run
```
	python <npdquast_installation_dir>/npd_quast.py total <work_folder>
```
