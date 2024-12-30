library(ctrialsgov)

# load all the data from clinicaltrials.gov into a cache (ca. 500MB)
# this takes some time 
ctgov_load_cache()

# query all phase 3 randomized studies with parallel assignment from the cache
data_load <- ctgov_query(
    phase = c("Phase 3"),
    allocation = "Randomized",
    intervention_model = "Parallel Assignment",
    )

# only select studies with the completion_date between 2013 and 2024
data <- data_load[data_load$completion_date >= "2004-01-01" & data_load$completion_date <= "2024-12-31", ]

# select only completed studies
data <- data[data$rec_status == "Completed", ]

# remove studies with unknown number of participants
data <- data[!is.na(data$enrollment), ]

# get the total number of studies
nrow(data) 

# get the median number of participants
median(data$enrollment) 
