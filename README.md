# ClassiCOL version 1.0.0
Code and usersguide

steps and explanations below will be addapted shortly

1. Download the zip folder named 'ClassiCol_version_1_0_0'
2. Unzip the file in a location of your chosing
3. Open prompt and change you r workdirectory to the location of the Classicol_1_0_0 python script
4. use the following commands to start the algorithm:
> python Classicol_version_1_0_0.py -d path_to_the_script\Classicol_version_1_0_0
> -l Demo (in case you want to test the algorithm) or the folder location containing your personal Mascot *.csv or maxquant *.txt output files
> -s MASCOT or MaxQuant (specify the search engine used)
> -t (optional) you can restrict the taxonomy by specifying it e.g. Pecora or multiple taxonomy Pecora|Primates
> -m specify the fixed modification used during protein extraction e.g. C,45.987721 or multiple with C,45.98|M,...
> -f (optional) location of the folder containing a custom database in fasta format
WARNING: The algorithm can use a substantial amount of the available CPU's and memmory. When not enough is free there is a chance the algorithm will go into error
