load_rdata = function() { 
  
  files <- list.files(here::here("data", "raw_data"), 
                      pattern = ".RData|Rdata", 
                      full.names = TRUE)

  load_path <- lapply(files, 
                      load, 
                      .GlobalEnv)  
  
  files_csv = list.files(here::here("data", "raw_data"), 
                         pattern = ".csv", 
                         full.names = T)
  
  load_csv = lapply(files_csv,
                    function(i){read.csv(i, header = TRUE)})
  names(load_csv) = c("covar_graham", "covar_jb")
  list2env(load_csv, envir = .GlobalEnv)
  
  }