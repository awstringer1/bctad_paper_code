# BCT AD Paper Code
Code for the BCT AD paper. All required packages and data are downloaded within the scripts. 

Currently, user has to change the paths in the scripts to match those on their machine. The following steps can be followed:

1. Unzip/clone this repository at your preferred location on your machine,
2. In scripts `01-relative-accuracy.R`, `02-hodges.R` and `03-hodges-heavytailed.R`, change the `global_path` variable to the directory where you unzipped/cloned this repository.
3. Run those scripts. The results of Sections 2 and 3 (Tables 1 and 2 and Figures 1 and 2), respectively, will appear in the `figures` folder. Script `02-hodges.R` takes some time to run, because of the model-fitting timing comparisons; on my hardware it took about 30 minutes.
