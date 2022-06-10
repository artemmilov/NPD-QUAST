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
