library(ctrialsgov)

# loads all the data from clinicaltrials.gov into a cache (ca. 500MB)
# this takes time 
ctgov_load_cache()

# query all phase 3 studies from the cached data
data_load <- ctgov_query(
    phase = "Phase 3",
    allocation = "Randomized",
    intervention_model = "Parallel Assignment",
    primary_purpose = "Treatment",
    )


# only select studies with the completion_date year in 2023
data <- data_load[data_load$completion_date >= "2013-01-01" & data_load$completion_date <= "2024-12-31", ]

# remove studies with unknown number of participants
data <- data[!is.na(data$enrollment), ]

nrow(data) # number of studies
median(data$enrollment) # median number of participants (if known)

