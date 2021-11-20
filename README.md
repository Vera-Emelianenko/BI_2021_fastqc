# BI_2021_fastqc

This is python project for python course in [Bioinnformatics Institute](https://bioinf.me/en). The fastqc.py is a script for sequencing reads quality control. It provides the similar metrics and graphs as the widely used FastQC tool from Babraham Institute. 

## Contributors



## Installation

The installation process does not rely on distinct OS features and can be applied for Windows, Linux and Mac systems. The script is supposed to be run from bash-like command line. To use it on Windows, please install Windows Subsystem for Linux (WSL) (this [user guide](https://www.windowscentral.com/install-windows-subsystem-linux-windows-10) can help you). Then you can run the commands in the WSL command line. The script is written and tested in Python3.8, and the installation process assumes you have python3 on your mashine. You can check the current version of python using `python --version`. To install or update Python, visit https://www.python.org/downloads/.

#### Downloading
You can download the archive with fastqc.py and test data using wget:

```
wget https://github.com/Vera-Emelianenko/BI_2021_fastqc/archive/refs/heads/main.zip
unzip main.zip
cd BI_2021_fastqc-main
```
Alternatively, the archive with the content of the project can be downloaded via visual interface of GitHub. Go to *main* branch, choose green button "Code" on the upper right and download zip. Zip archive might be unpacked with any tool you have or with unzip (`unzip BI_2021_fastqc.zip`)

You can also download the project via git clone: 
```
git clone https://github.com/Vera-Emelianenko/BI_2021_fastqc.git
````

### Dependencies installation

Since fastqc.py requires some particular versions of libraries, listed in requirement.py, the easiest way to run the script is to create virtual environment first and to install all necessary dependencies there. 

Create virtual environment:
```
python3.8 -m pip install --upgrade pip setuptools virtualenv
python3.8 -m virtualenv pseudofastqc # you can give any name you like 
source pseudofastqc/bin/activate
```
Note: If you are using not the bash shell, you may want to use some other file listed in pseudofastqc/bin, for example to activate the environment in fish, run 
`source fastqc/bin/activate.fish`

When environment is activated, you will see its name in brachets somewhere in your command line, e.g.
```
(pseudofastqc) vera@LAPTOP-EAJJI10O /m/c/U/e/D/b/BI_2021_fastqc-main>  
```

Install all dependencies: 
```
python -m pip install -r requirements.txt
```

You should see the following message: 
```
Successfully installed biopython-1.79 cycler-0.11.0 kiwisolver-1.3.2 matplotlib-3.4.3 numpy-1.21.4 pandas-1.3.4 pillow-8.4.0 pyparsing-3.0.6 python-dateutil-2.8.2 pytz-2021.3 scipy-1.7.2 seaborn-0.11.2 six-1.16.0
```
Alternatively, you can install all the dependencies in requirements.txt manually. To install any module, run the following:

`pip3 install module`

Note: in case you decided to run it in Python2, we can not quarantee proper script work. In this case the modules have to be installed with `pip install module`

### Verifying the installation

To verify that the script works, you can run it on the test data: 

```
python fastqc.py -i test_data/test.fastq -o ./
```
If you see no error messages, you are all set! You can remove test results with `rm ./test*`. 

## Usage

To run fastqc from the command line, type

``` bash
python fastqc.py -i input_file_dir/input_file.fq -o output_file_dir
```

## Basic options

`-i <path_to_file>` (or `--input <path_to_file>`)
    Reads data in fastq format.
    
`-o <output_dir>` (or `--outdir <output_dir>`)
    Specify the output directory. The current directory is default `"./"`.

## Output files 


## Operating system and software versions

