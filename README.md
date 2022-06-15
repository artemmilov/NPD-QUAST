# NPD-Quast
NPD-Quast is designed for analyzing natural peptide identificator\`s reports of the provided data and comparing them to each other.

* Copability: NPD-Quast is Python 3.9+ compliant
* Prerequisites: Apart from a working version of Python, the following non-standard libraries need to be installed: conda, numpy and matplotlib. Also you need to install RDKit if you want to run NPD-Quast whithout `npd-quast-env` activated.
* Disclaimer: NPD-Quast should work on all platforms, but was only tested in a Linux (Ubuntu 20.04) system.

## Installing NPD-Quast
Just clone this repository on your local machine and fill `npd_quast.ini` file. Note that you need to write location just that tools, which you want to run.

## Running NPD-Quast
To run NPD-Quast you firstly need to construct work folder which structure is written below. Put all raw data in "challenges/" dir and all tool answers files in "reports/" dir. Also, you need all dependencies to be installed. We recommend you activate a "npd-quast-env.yml" environment which is located in "envs/" dir. So, just run this command:

```conda activate npd-quast-env```

After that, if you have putted tool answers in "reports/" folder than you need to make a full reports from them. To do this you have to run such command:

```python npd_quast.py compile_reports your_work_folder_name```

So now you have a site with pages corresponding to your reports. If you want to add new report from unsupported tool. You have to make tool answers by yourself and make new dir with them in "reports/" dir. Next, run command "compile_reports" again.

Also, if you want to add new report from supported tool. Just tap this command.

```python npd_quast.py run_n_report tool_name report_folder your_work_folder_name```

## Work folder structure
The work folder consist of directory with challenges,
directory with report and file with true answers.

True answers are written as strings <i>"challenge_name$spectra_name \t answer_inchi"</i>,
where <i>answer_inchi</i> is the first 14 symbols of inchi-key of molecule corresponding to
the spectrum <i>spectra_name</i> from the challenge <i>challenge_name</i>.

This is followed by a directory with challenges.
Each challenge consist of database and spectra list. Database should be written
in <i>csv</i> format and all spectra in <i>mgf</i>. Also each spectrum file
contain only one spectrum. The tool under test makes identifications between
spectra and molecules from database, corresponding to this challenge.

The <i>"reports"</i> directory contains pages of this report and
answers of all tested tools. All identifications, which tool with
<i>tool_name</i> has conducted are located in file
<i>"./reports/tool_name/tool_answers.txt"</i>.
They are written as strings <i>"challenge_name$spectra_name \t answer_inchi"</i>,
where <i>answer_inchi</i> is the first 14 symbols of inchi-key of molecule
which this tool has matched to the spectrum <i>spectra_name</i> of
<i>challenge_name</i> challenge. <i>"scan"</i> is a number,
which corresponds to the confidence of match.
Tested tool is more confident if <i>"scan"</i> is less.

Based on the tool answers and true answers NPD-Quast makes a report on the quality of work
each tool on all data and also comparing them to each other.
