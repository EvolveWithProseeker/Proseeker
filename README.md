# Proseeker

In silico directed evolution

Proseeker is supplied as a standalone Python script (Proseeker.py) for either interative use or pipeline integration. 

If running in a pipeline the commenting on the below lines (38 & 40 in Proseeker.py as written) shouled be uncommented and the path to the working directory supplied as a commandline argument:

  user_input = input("Enter the path to your Proseeker working directory (i.e C:\Proseeker): ")

  # user_input = int(sys.argv[1])

The Proseeker working directory should contain the files ranking.csv, jobstart.csv and whatever you have named your Bres file directory (i.e. PPRONLY).

The ranking file contains values drawn from AAindex as stated in the paper.

The jobstart.csv file (example included) is laid out as below in two lines with headers as shown:

  g,k,MEG,block,g1,pchoice,A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,truncres,cores,sites,bres
  20,100,2,MMMMFFFFFFFFMMMMMMMMMMMMMFFFFFFFFFFFFFMM,SELVNEIVKQLAEVAKEATDKELVIYIVKILAELAKQSTD,1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,10,4,2,BRES

g = number of generations and should be any positive integer greater than 1 (20 is the default).
k = the number of variants to be amde per generation. This should be a positive integer greater than 10 (100 is the default and is recommended).
MEG = the number of mutational events to occur per generation (this differs from EMEG as stated in the paper).
block = the masking scheme to be used. M=mutable and F = fixed. This is not option so if no masking is desired simply include M at all positions. This string must be the same length as the progenitor protein.
g1 = the progenitor protein. This should contain only a string consisting of letters (case is irrelevant) corresponding to one of the 20 natural amino acids. STOP codon is not included and Proseeker has no ability to introduce STOP codons.
pchoice = the choice to use a specified probability set or not. 1 = yes and 0 = no. 
A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V = the individual probability of selection of each AA. These should add up to 1 and as such 0.05 represent an even probability of selection. If phoice = 1 then the values listed here will be used. 
truncres = Proportion of results to provide to the user. The default is 10 which will results in the best k/10 ressults per generation being provided to the user.
cores = specification of number cores. This is present to support the parallelized version. This is set to 4 default however for non-cluster use it may be ignored.
sites = The number of sites per protein within the assessment library that are to be used as comparators
BRES = the name of the directroy in which your BRES files are stored (i.e. if the BRES files are in C:\path\to\Proseeker\BRES then you would specify just BRES.

Once run Proseeker will create a time-stamped directory for the run within the working directory. This directroy will be empty at first however it will save each generations best results in the created file for later harvesting and examination. As such multiple runs can be performed simultaneously on a desktop enironment allowing triplicate runs for verification and assessment of consistency.

# BresMaker

Making BRES files is acomplished using the BresMaker.py script.

This script may be run in a pipeline or in an interactive environment by changing commenting on the below lines to produce ether a commandline or interactive version (24:34 in BresMaker.py as written):

  #user_input = str(sys.argv[1])
  #ranking = str(sys.argv[2])
  #working = str(sys.argv[3])
  #iterations = int(sys.argv[4])
  #trys = int(sys.argv[5])

  user_input = "D:/Proseeker/exampledeets.csv"
  ranking = "D:/Proseeker/ranking.csv"
  working = "D:/Proseeker"
  iterations = 1000000
  trys = 1000

The details file (exampledeets.csv) loaed as 'user_input' consists of a .csv file with two columns. column one consists of the sequence of the protein from which sites are to be harvested as below:

  LTPEQVVAIASNIGGKQALETVQALLPVLCQAHG,12:13:14

Column 1 consists of one string per row consisting of letters (case is irrelevant) corresponding to one of the 20 natural amino acids. STOP codon is not included and Proseeker has no ability to introduce STOP codons.

Column 2 consists of positive integers (you could use negative integers as Python will take -1 for example to refer to the last character however it is easy to get negative indexing wrong). If only 1 site is to be taken then 1 integer is provided. If more than 1 site is required then 1 integer per site is required separated by : with no spaces.

The ranking file (ranking.csv) loaded as 'ranking' contains values drawn from AAindex as stated in the paper.

Working is the directory in which the BresFiles are to be created. It is recommended that a specific directory be created for this. 

Iterations and Trys refer to the maximum number of k-means clustering interations to be performed and the number of time the k-means algorithm will be run with different centroid seeds respectively. 
