check_var1_vs_multiple_var2 <- function(var1, var2_vector, df) {
  # Check if var1 exists in the dataframe
  if (!(var1 %in% names(df))) {
    stop(paste(var1, "column is missing from the dataset."))
  }
  
  # Initialize a list to store results
  results <- list()
  
  for (var2 in var2_vector) {
    # Check if var2 exists in the dataframe
    if (!(var2 %in% names(df))) {
      warning(paste(var2, "column is missing from the dataset. Skipping."))
      next
    }
    
    # Check if there are non-NA values in both columns
    valid_data <- sum(!is.na(df[[var1]]) & !is.na(df[[var2]]))
    
    if (valid_data == 0) {
      message(paste("There are no valid (non-NA) pairs of data for", var1, "and", var2, "."))
    } else {
      message(paste("Found", valid_data, "valid pairs of data for", var1, "and", var2, "."))
      results[[var2]] <- valid_data
    }
  }
  
  if (length(results) == 0) {
    message("No valid pairs found for any variables.")
    return(NULL)
  } else {
    return(results)
  }
}

get_response_description <- function(a, response_variable) {
  # Check if the required data frame exists in a
  if (!("mesi_response_variable_abbreviations" %in% names(a))) {
    stop("The 'mesi_response_variable_abbreviations' data frame is not present in the provided object.")
  }
  
  # Get the description
  description <- a$mesi_response_variable_abbreviations$description[
    a$mesi_response_variable_abbreviations$response == response_variable
  ]
  
  # Check if a description was found
  if (length(description) == 0) {
    warning(paste("No description found for response variable:", response_variable))
    return(NULL)
  } else if (length(description) > 1) {
    warning(paste("Multiple descriptions found for response variable:", response_variable))
    return(description)
  } else {
    return(description)
  }
}

mesi_import <- function() {
  mesi_path = "/home/josimms/Documents/Austria/MESI_data/data"
  
  a <- lapply(list.files(mesi_path, pattern = "mesi", full.names = TRUE), data.table::fread)
  names(a) <- gsub(".csv", "", list.files(mesi_path, pattern = "mesi"))
  
  # Key!
  get_response_description(a, "npp")
  # Posibilities
  posibilities = unique(a$mesi_response_variable_abbreviations$description)

  # Transform the data so the response variables are columns
  df <- a$mesi_main %>%
    group_by(id, response) %>%
    pivot_wider(names_from = response, values_from = x_c, values_fn = mean)
  
  # TODO: none of the variables have any matches, so need to fix this!
  check_var1_vs_multiple_var2("npp", posibilities, df)
  
  ###
  # Boreal
  ###
  boreal_main_data = a$mesi_main[a$mesi_main$vegetation_type == "boreal_forest",]
  
  # TODO: need to transform
  
  # Posibilities
  boreal_posibilities = unique(boreal_main_data$response)
  check_var1_vs_multiple_var2("jmax", boreal_posibilities, boreal_main_data)
  
  
  return(boreal_main_data)
}

try_data <- function() {
  # Read in the data
  direct = "/home/josimms/Documents/Austria/TRY_data/"
  data_try <- data.table::fread(paste0(direct, "35495.txt"))
  
  # Seperate for pine
  data_try_pine <- data_try[data_try$SpeciesName == "Pinus sylvestris",]
  
  # Make the column of variables into lots of columns
  out_data_try_pine_pivot <- data_try_pine %>%
    filter(TraitName != "") %>% # TODO: check why there are empty TraitName rows
    # Group by ObservationID and TraitName, and summarize StdValue
    # TODO: why are the multiple values here?
    group_by(ObservationID, TraitName) %>%
    summarize(StdValue = mean(StdValue, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = TraitName, values_from = StdValue)
  
  # TODO: check that this is working properly
  
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


