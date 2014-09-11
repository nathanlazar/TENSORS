# nathan dot lazar at gmail dot com

# Functions used in manipulating tensor data

library(dplyr)

# Function to use dplyr to average responses for each dose
avg_resp <- function(resp, cl, cpd) {
tmp <- resp %.%
  filter(cellline==cl) %.%
  filter(compound==cpd) %.%
  select(od0.1:od9.3)
  if(nrow(tmp)>0) {
    tmp <- tmp %.%
    summarise_each(funs(mean)) %.%
    mutate(od0 = mean(od0.1, od0.2, od0.3),
           od1 = mean(od1.1, od1.2, od1.3),
           od2 = mean(od2.1, od2.2, od2.3),
           od3 = mean(od3.1, od3.2, od3.3),
           od4 = mean(od4.1, od4.2, od4.3),
           od5 = mean(od5.1, od5.2, od5.3),
           od6 = mean(od6.1, od6.2, od6.3),
           od7 = mean(od7.1, od7.2, od7.3),
           od8 = mean(od8.1, od8.2, od8.3),
           od9 = mean(od9.1, od9.2, od9.3)) %.%
    select(od0:od9)
    return(as.numeric(tmp))
  } else {
    return(rep(NA, 10))
  }
}

