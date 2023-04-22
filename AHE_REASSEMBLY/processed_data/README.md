pre-print available: https://doi.org/10.1101/2022.10.21.513146

## Jacob S. Berv<sup>1,2,3,*</sup>, Sonal Singhal<sup>4</sup>, Daniel J. Field<sup>5,6</sup>, Nathanael Walker-Hale<sup>7</sup>, Sean W. McHugh<sup>8</sup>, J. Ryan Shipley<sup>9</sup>, Eliot T. Miller<sup>10</sup>, Rebecca T. Kimball<sup>11</sup>, Edward L. Braun<sup>11</sup>, Alex Dornburg<sup>12</sup>, C. Tomomi Parins-Fukuchi<sup>13</sup>, Richard O. Prum<sup>14,15</sup>, Benjamin M. Winger<sup>1,3</sup>, Matt Friedman<sup>2, 16</sup>, Stephen A. Smith<sup>1</sup>

Corresponding Author: *Jacob S. Berv, jberv@umich.edu, jacob.berv@gmail.com
Relevant co-author: Sonal Singhal, sonal.singhal1@gmail.com

Author affiliations:

1.	Department of Ecology and Evolutionary Biology, 1105 North University Avenue, Biological Sciences Building, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA

2.	University of Michigan Museum of Paleontology, 1105 North University Avenue, Biological Sciences Building, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA

3.	University of Michigan Museum of Zoology, 1105 North University Avenue, Biological Sciences Building, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA

4.	Department of Biology, California State University, Dominguez Hills, Carson, California 90747, USA

5.	Department of Earth Sciences, Downing Street, University of Cambridge, Cambridge CB2 3EQ, UK

6.	Museum of Zoology, Downing Street, University of Cambridge, Cambridge CB2 3EJ, UK

7.	Department of Plant Sciences, Downing Street, University of Cambridge, Cambridge, CB2 3EA, UK

8.	Department of Evolution, Ecology, and Population Biology, Washington University in St Louis, St Louis, Missouri, USA

9.	Department of Forest Dynamics Swiss Federal Institute for Forest, Snow, and Landscape Research WSL Zürcherstrasse 111 8903 Birmensdorf, Switzerland

10.	Macaulay Library, Cornell Lab of Ornithology, Ithaca, New York, 14850, USA

11.	Department of Biology, University of Florida, Gainesville, Florida 32611, USA

12.	Department of Bioinformatics and Genomics, University of North Carolina at Charlotte, Charlotte, North Carolina, USA

13.	Department of Ecology and Evolutionary Biology, University of Toronto, Toronto, Ontario, Canada, M5S 3B2

14.	Department of Ecology and Evolutionary Biology, Yale University, New Haven, Connecticut, 06520, USA

15.	Peabody Museum of Natural History, Yale University, New Haven, Connecticut, 06520, USA

16.	Department of Earth and Environmental Sciences, 1100 North University Avenue, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA


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
	- there were two runs to re-do some of the failed samples from run1 and to add additional species
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