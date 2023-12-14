vcsn_info <- function(VCSN.directory) {

  # Find the year range for the VCSN data
  #    - zip file names follow the format:  vcsn_1960.zip
  zip.files <- list.files(VCSN.directory, pattern = ".zip")
  vcsn.years <- as.numeric(substr(zip.files, start = 6, stop = 9))

  info <- NULL
  info$VCSN.directory <- VCSN.directory
  info$first.zip.file <- zip.files[1]
  info$min.year <- min(vcsn.years, na.rm = TRUE)
  info$max.year <- max(vcsn.years, na.rm = TRUE)

  return(info)

}


vcsn_agent_locations <- function(VCSN.directory) {
  #
  # Find the virtual lat/lon coordinates and associated Agent number
  #    -  it is assumed they are the same for all days and years
  #    -  so they are pulled for the first day in 1960
  #       (i.e. first file within the earliest zip file )
  #
  # Arguments:
  #
  # VCSN.directory: character string for the VCSN directory location
  #                  e.g.  "Q:/CLIMATE/vcsn_data/"
  #
  # Value:
  #
  # A sf object with four columns: Agent (numeric), vcsn.lon, vcsn.lat geometry
  #
  # where the CRS is NZGD49 (EPSG: 27258)??
  #
  # Assume CRS is 4326

  vcsn.info <- vcsn_info(VCSN.directory)

  earliest.zip.file.name <- paste0(VCSN.directory, "vcsn_", vcsn.info$min.year, ".zip")
  fnames <- as.character(unzip(earliest.zip.file.name, list = TRUE)$Name)

  first.daily.file <- read_csv(unzip(zipfile = earliest.zip.file.name, files = fnames[1]),
                               col_types = cols())

  # 11491 virtual lat/lon locations
  #    - virtual lat/lon are to three decimal places
  vcsn.locations <- first.daily.file %>%
    dplyr::select(Agent, Longt, Lat) %>%
    rename(vcsn.lat = Lat, vcsn.lon = Longt) %>%
    distinct() %>%
    st_as_sf(coords = c("vcsn.lon", "vcsn.lat"), remove = FALSE, crs = 4326)

  return(vcsn.locations)

}


vcsn_append_nearest_vcsn <- function(input.df, 
                                     VCSN.directory, 
                                     input.df.EPSG = 4326) {

  # Find the "Agent" (i.e virtual climate station) closest to
  # the points in the data frame input.df and append these to the
  # data frame
  #
  # Arguments:
  #
  # input.df: data frame with columns lon, lat
  #
  # VCSN.directory: character string with directory name containing vcsn zip files
  #                    eg. "Q:/CLIMATE/vcsn_data/"
  #
  # input.df.EPSG:  EPSG number for the lon, lat coordinate reference system with a
  #                 default of 4326 (i.e. WGS 84)
  #
  #
  # Value: original input data frame with additional columns:
  #
  #
  #
  

  if(!("lon" %in% names(input.df))) stop("You need a lon column for the input data frame")
  if(!("lat" %in% names(input.df))) stop("You need a lat column for the input data frame")

  if(missing(input.df.EPSG)) warning("The assumed EPSG is 4326 for your input data frame\n")

  # Find all the virtual climate station locations
  vcsn.locations <- vcsn_agent_locations(VCSN.directory)

  # This needs a fix-up, a transform of the coordinate reference
  # system, instead of simply labelling with
  input.sf <- input.df %>%
    st_as_sf(coords = c("lon", "lat"), remove = FALSE, crs = input.df.EPSG)
    # st_transform(crs = st_crs(vcsn.locations)) %>%
    # st_transform(crs = "+proj=longlat")

  nearest.Agent.index <- input.sf %>%
    st_set_crs(st_crs(vcsn.locations)) %>%
    st_nearest_feature(vcsn.locations)

  vcsn.locations.nearest <- vcsn.locations %>%
    slice(nearest.Agent.index)

  output.df <- input.df

  output.df$VCSN.Agent <-vcsn.locations.nearest$Agent

  nearest.positions <- vcsn.locations.nearest %>%
    st_coordinates() %>%
    as.data.frame() %>%
    rename(VCSN.lon = X, VCSN.lat = Y) %>%
    relocate(VCSN.lon, .after = VCSN.lat)

  output.df <- bind_cols(output.df, nearest.positions)

  return(output.df)
}


vcsn_append_climate_information <- function(input.df,
                                            VCSN.directory,                                     
                                            years) {

  input.sf.plus.vcsn <- vcsn_append_nearest_vcsn(input.df, VCSN.directory)

  # Loop over years first, so only pulling out each zip file once

  all.ev.data <- NULL
  tempdir <- tempdir()

  for (zip.year in years) {

    message("Processing VCSN data for year: ", zip.year, "\n")

    # clear temporary directory first
    unlink(paste0(tempdir, "\\", "*vcsn.dat"))

    # unzip into temporary directory
    zip.file.name <- paste0(VCSN.directory, "vcsn_", zip.year, ".zip")
    unzip(zip.file.name, exdir = tempdir)

    # list all the unzipped vcsn data files
    daily.file.names <- list.files(tempdir, pattern = glob2rx("*vcsn.dat"), full.names = TRUE)

    for (this.day.file in daily.file.names) {

      message("Processing day file ", this.day.file, " for VCSN data")

      daily.file <- read_csv(this.day.file, col_types = cols())

      # data files at the very end may be empty
      if(length(daily.file) > 0) {
        for (i in 1:nrow(input.sf.plus.vcsn)) {
          this.Agent <- input.sf.plus.vcsn$VCSN.Agent[i]
          this.ev <- filter(daily.file, Agent == this.Agent)
          this.agent.ev <- bind_cols(input.sf.plus.vcsn[i, ], this.ev)
          all.ev.data <- bind_rows(all.ev.data, this.agent.ev)

        }
      }

    }

  }

  # Already in data (VCSN.Agent, VCSN.lat, VCSN.lon) so remove
  all.ev.data <- all.ev.data %>%
    dplyr::select(-Agent, -Lat, -Longt)

  # Convert from character to Date class
  all.ev.data$Date = as.Date(all.ev.data$Date, format = '%d/%m/%Y')

  return(all.ev.data)

}



vcsn_combine_daily_data<- function(VCSN.directory, years) {
  
  # Loop over years first, so only pulling out each zip file once
  
  combine.vcsn.data <- NULL
  tempdir <- tempdir()
  
  for (zip.year in years) {
    
    message("Processing VCSN data for year: ", zip.year, "\n")
    
    # clear temporary directory first
    unlink(paste0(tempdir, "\\", "*vcsn.dat"))
    
    # unzip into temporary directory
    zip.file.name <- paste0(VCSN.directory, "vcsn_", zip.year, ".zip")
    unzip(zip.file.name, exdir = tempdir)
    
    # list all the unzipped vcsn data files
    daily.file.names <- list.files(tempdir, pattern = glob2rx("*vcsn.dat"), full.names = TRUE)
    
    for (this.day.file in daily.file.names) {
      message("Processing day file ", this.day.file, " for VCSN data")
      daily.file <- read_csv(this.day.file, col_types = cols())
      # data files at the very end may be empty
      if(length(daily.file) > 0) {      
          combine.vcsn.data <- bind_rows(combine.vcsn.data, daily.file)
      }
    }
    
  }
  
  # Convert from character to Date class
  combine.vcsn.data$Date = as.Date(combine.vcsn.data$Date, format = '%d/%m/%Y')
  
  return(combine.vcsn.data)  
  
}


