mesi_import <- function() {
  mesi_path = "/home/josimms/Documents/Austria/MESI_data/data"
  
  a <- lapply(list.files(mesi_path, pattern = "mesi", full.names = TRUE), data.table::fread)
  names(a) <- gsub(".csv", "", list.files(mesi_path, pattern = "mesi"))
  
  boreal_main_data = a$mesi_main[a$mesi_main$vegetation_type == "boreal_forest",]
  
  return(boreal_main_data)
}