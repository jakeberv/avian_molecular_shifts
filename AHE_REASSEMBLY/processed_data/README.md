
## Information for Processed Data 
---

- structure of this folder

	```
	- processed data
		- run1
			- PRG
			- samples1.csv
			- variants
		- run2
			- PRG
			- samples2.csv
			- variants
	```

- run1 _versus_ run2
	- run1 consists of 196 samples
	- run2 consists of 22 samples
	- For Prum et al 2015, there were two sequencing runs to re-do some of the failed samples from run1 and to add additional species
- information about these files
	- `samples1.csv` and `samples2.csv`
		- contains information about the species names, their sample name, and their adaptor sequences
	- PRG
		- PRG: pseudo-reference genome
		- these files should be one per species
		- they consist of annotated loci (in coding orientation)
		- the number of loci in each file will be slightly variable depending on how many loci were assembled
	- variants
		- these are variant call files (VCFs)
		- per species, there are three files
			- files include variable & invariable sites
			- all have sites filtered to only include those with quality >20
			- these were filtered further by coverage, only including variants with >2, >5, >10x coverage