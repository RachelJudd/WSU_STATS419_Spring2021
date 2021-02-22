###  ASSIGNMENT (1) ###
You will need to read enough of the article to understand what those datasets are. [Wells-23.pdf]

You will need to figure out how to read in these two "datasets" found in this PDF.  "Read in" is a generic term, it does not demand a function, this is a "git r done" tasks, maybe not very glamorous.  Copy/Paste, whatever works.

Table 1 & 2 has information on Wells 1-23 and auxilliary data (mean, min/max, etc.)  Your dataframes should only contain the information for the 23 wells.  The first dataframe 'owell.metals' [Table 1] should initally have 15 columns (and 23 rows).  The second dataframe 'owell.chemistry' [Table 2] should initially have 20 columns (and 23 rows).

I am providing an auxiliary file 'wells-location.txt' where I manually placed the 23 wells based on latitude, longitude, and altitude.  You will need to append this data to your two "dataframes".

You will need to manually append a new column to each data frame called "fault" (see Fig 1, pg4 of PDF [770]) and estimate in "meters" a value for each well as a nearest distance to a fault.

You will need to manually append a new column to each data frame called "geology" (see Fig 1, pg4 of PDF [770]) and assign a comma-separated string of the first word of relevant rock formations near the wells.  For example, Well-23 has a pink and salmon color, so its geology value would be:  "Alkaline,Granite"

Create a list with these two dataframes:  "wells$saudi$chemistry" and "wells$saudi$metals" ... save this list as RDS named "wells.rds"

Create a notebook documenting and describing the datasets, calling it 'oWell-1.Rmd' ...

This is o-well assignment #1.

###  ASSIGNMENT (2) ###

- Let's build one datafile "saudi.txt" that has one big data frame.


Perform some basic analyses and graphics summarizing the data in a notebook 'oWell-2.Rmd' ....  think about the data content available and any curious questions you may have.

Descriptive statistics is generally a univariate procedure (recall your Stats 360?) that can be applied to multiple variables. 

This is o-well assignment #2. 


