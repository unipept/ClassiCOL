# ClassiCOL version 1.0.0
## Code and User's Guide

The steps and explanations below will be adapted shortly.

### Installation

1. Download the code in this repository.
2. Navigate to the download location using the command prompt.
3. Install the required packages using `pip install -r requirements.txt`.

### Usage
Use the following command to start the algorithm with the demo data:
```sh
$ python Classicol.py -d path_to_the_script -l Demo -s MASCOT
```

You can use the arguments as follows:
  - `-l Demo` (in case you want to test the algorithm) or the folder location containing your personal Mascot \*.csv or MaxQuant \*.txt output files
  - `-s MASCOT` or `MaxQuant` (specify the search engine used)
  - `-t` (optional) you can restrict the taxonomy by specifying it, e.g., Pecora or multiple taxonomy Pecora|Primates
  - `-m` specify the fixed modification used during protein extraction, e.g., C,45.987721 or multiple with C,45.98|M,...
  - `-f` (optional) location of the folder containing a custom database in fasta format

**WARNING:** The algorithm can use a substantial amount of the available CPU and memory. When not enough is free, there is a chance the algorithm will go into error.
