mesi_import <- function() {
  mesi_path = "/home/josimms/Documents/Austria/MESI_data/data"
  
  a <- lapply(list.files(mesi_path, pattern = "mesi", full.names = TRUE), data.table::fread)
  names(a) <- gsub(".csv", "", list.files(mesi_path, pattern = "mesi"))
  # Key!
  a$mesi_response_variable_abbreviations$description[a$mesi_response_variable_abbreviations$response == "anpp"]
  # Posibilities
  unique(a$mesi_response_variable_abbreviations$description)
  
  ###
  # Boreal
  ###
  boreal_main_data = a$mesi_main[a$mesi_main$vegetation_type == "boreal_forest",]
  
  # TODO: go through this data when I have the life history model working!
  unique(boreal_main_data$response)
  
  ###
  # Jmax
  ###
  main_data_jmax <- a$mesi_main[a$mesi_main$response %in% c("jmax", "vcmax", "leaf_nitrogen_concentration_by_mass"),]
  main_data_jmax
  
  main_data_jmax
  
  a$mesi_main[a$mesi_main$response %in% c("leaf_n_area", "jmax"),]
  
  plot(a$mesi_main$x_c[a$mesi_main$response == "leaf_n_area"])
  
  df = pivot_wider(a$mesi_main, 
                   names_from = response, 
                   values_from = x_c,
                   values_fn = mean)
  
  a$mesi_main$x_units[a$mesi_main$response == "leaf_n_mass"]
  
  return(boreal_main_data)
}

n_content_data <- function() {
  direct = "/home/josimms/Documents/Austria/extra_data/"
  
  # https://zenodo.org/records/11066705
  data <- readxl::read_excel(paste0(direct, "Muhonen et al article metadata file_ UEF Peltola Heli 25042024.xlsx"), sheet = 1)
  data_scots_pine = data[data$`Main tree species` == "Scots pine",]
  
  par(mfrow = c(1, 2))
  plot(data_scots_pine$`Fertilization dose (kg N/ha)`, data_scots_pine$`N needles (%) in 2021`, xlab = "Fertilization Dose", ylab = "N Needles %")
  plot(data_scots_pine$`N humus (%) in 2022`, data_scots_pine$`N needles (%) in 2021`, pch = "x", col = "blue", xlab = "N humus (%) in 2022", ylab = "N needles (%) in 2021")
  
  ### TODO: Missing / would be useful
  # Nitrogen per area to nitrogen per mass
  
  data_collation <- data.table::data.table(
    site = c("Hyytiälä",
             "France", 
             "France",
             "Mekrijärvi Research Station ",
             "Mekrijärvi Research Station ",
             "Mekrijärvi Research Station ",
             "Mekrijärvi Research Station "),
    reference = c("Korhonen, 2007", 
                  "Warren, 2002",
                  "Warren, 2002",
                  "Kellomäki, 1999",
                  "Kellomäki, 1999",
                  "Kellomäki, 1999",
                  "Kellomäki, 1999"),
    n_content = c(1.21, NA, NA, NA, NA, NA, NA),
    a_jmax = c(NA, NA, NA, 32.45, 23.24, 21.46, 30.28),
    a_vcmax = c(NA, 11.319, 4.74, 15.36, 14.91, 12.06, 12.74),
    details = c(NA,
                "Current year, Jmax Vcmax closely and linearly related",
                "1 year old, Jmax Vcmax closely and linearly related",
                "Elevated C, nitrogen leaf per area",
                "Elevated C and T, nitrogen leaf per area",
                "Elevated T, nitrogen leaf per area",
                "Control, nitrogen leaf per area")
  )
  
}


dFdIb_plot <- function() {
  dFdIb <- function(Ib, alpha, a_jmax, N) {
    return((alpha * a_jmax * N) / Ib^2)
  }
  
  Ib_range = seq(0.75, 2, length.out = 50)
  plot(Ib_range, dFdIb(Ib_range, 0.2, 800, 0.02), col = "red", 
       type = "l", xlab = "Range of Ib values", ylab = "Marginal cost on F of Ib",
       ylim = c(0, max(dFdIb(Ib_range, 0.2, 800, 0.02))))
  lines(Ib_range, dFdIb(Ib_range, 0.2, 800, 0.01), col = "orange")
  lines(Ib_range, dFdIb(Ib_range, 0.2, 800, 0.005), col = "yellow")
  legend("topright", c("Leaf Nitrogen = 2%", "Leaf Nitrogen = 1%", "Leaf Nitrogen = 0.5%"), 
         col = c("red", "orange", "yellow"), bty = "n", lty = 1)
}


