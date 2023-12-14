
iso_master_show_sites <- function(path.file.location, data.type = "") {

  if (!(data.type %in% c("Daily", "Weekly", "Monthly"))) {
    stop('The argument data.type must be one of: "Daily", "Weekly", "Monthly"')
  }
  
  master.data <- readxl::read_xlsx(path.file.location, 
                                   sheet = data.type, 
                                   guess_max = 2500)  

  gps.categories <- dplyr::distinct(master.data,  Group, Project, Site)
  print(gps.categories, n = 200)
}

iso_master_extract_data <- function(path.file.location,
                                    data.type = "",
                                    Site.contains = "",
                                    Date.start = NA,
                                    Date.end = NA, 
                                    quiet = FALSE) {

  # path.file.location <- "data/master precipitation isotope data/precip.iso.database.nz.14th.April.2023.xlsx"
  # data.type = "Monthly"
  # Site.contains <- "FREW"
  
  if (!(data.type %in% c("Daily", "Weekly", "Monthly"))) {
    stop('The argument data.type must be one of: "Daily", "Weekly", "Monthly"')
  }
  
  master.data <- readxl::read_xlsx(path.file.location, 
                                   sheet = data.type, 
                                   guess_max = 2500)
  
  master.data <- master.data |>
    dplyr::filter(grepl(pattern = paste0("^", Site.contains), x= Site))
  
  # Convert strings to date objects e.g. "28\2\2020" -> to date object
  if(!is.na(Date.start)) Date.start.obj <- lubridate::dmy(Date.start)
  if(!is.na(Date.end))   Date.end.obj   <- lubridate::dmy(Date.end)
  
  # Just the start date given
  if(!is.na(Date.start) & is.na(Date.end)) {
    master.data <- master.data |>
      dplyr::filter(Date >= Date.start.obj)
  }
  
  # Just the end date given
  if(is.na(Date.start) & !is.na(Date.end)) {
    master.data <- master.data |>
      dplyr::filter(Date <= Date.end.obj)
  }  
  
  # Start and end dates given
  if(!is.na(Date.start) & !is.na(Date.end)) {
  master.data <- master.data |>
    dplyr::filter(dplyr::between(Date, Date.start.obj, Date.end.obj))
  }
 
  if(!quiet) {
    gps.combinations <- dplyr::distinct(master.data, Group, Project, Site)
    message("")
    message("unique Group, Project, Site combinations for extracted data are:")
    message("")
    print(gps.combinations, n = 200)
  }
 
  return(master.data)

}




iso_Julian_day_fraction <- function(Date) {
  # 
  # Date: a vector of date objects
  #
  # Return values:  fraction of the year for the date
  #                 e.g. 1st Jan = 1/365 or 1/366 (depending on the year)
  #
  #
  # e.g  
  # 
  # x <- c( make_date(2020, 12, 1), make_date(2020, 12, 2)) 
  # 
  # iso_Julian_day_fraction(x)
  #
  # gives  0.9180328 0.9207650
  #
  
  Year <-     lubridate::year(Date)
  Year.day <- lubridate::yday(Date)
  Year.num.days <- yday(lubridate::make_date(year = Year, day = 31, month = 12))
  fraction.year <- Year.day/Year.num.days
  
  return(fraction.year)
  
}


iso_fit_sinusoidal_model <- function(the.data,
                                     site.column.name,
                                     date.column.name,
                                     quantity.column.name,
                                     lower.limits = c(Amplitude = 0.1, phase = -3.1416, offset = -10),
                                     upper.limits = c(Amplitude = 5.0, phase =  3.1416, offset = -3)) {
  # e.g.
  #
  # ddata.AK1 <- filter(ddata, Site == "AK1")
  #
  # modAK1 <- iso_fit_sinusoidal_model(ddata.AK1, "Site", "Date", "d18O")

  # Julian day. 1st Jan = 1,
  #            last day of year = 365 (non-leap year) or 366 (leap year)
  #
  # lower.limit, upper.limit values are what was used for d18O
  #


  isotope.data <- the.data |>
    dplyr::select(
      Site = !!site.column.name,
      Date = !!date.column.name,
      quantity = !!quantity.column.name
    ) |>
    mutate(Julian.Day.Fraction = iso_Julian_day_fraction(Date), .after = Date)

  # Nest the data by site. Each row is a single site (with the data)
  by_site <- isotope.data |>
    group_by(Site) |>
    tidyr::nest()

  # Model to be fitted to a single site
  site_model <- function(df) {
    form <- quantity ~ Amplitude * sin(2 * 3.141593 * Julian.Day.Fraction - phase) + offset
    pnms <- c("Amplitude", "phase", "offset")


    offset.start.value <- mean(isotope.data$quantity, na.rm = TRUE)

    # From help for nlrob:
    #
    # For (the default) method = "M", if the bounds are unspecified all parameters
    # are assumed to be unconstrained; also, for method "M", bounds can only
    # be used with the "port" algorithm
    #

    # For methods "CM" and "mtl", the bounds must additionally have an entry
    # named "sigma" as that is determined simultaneously in the same optimization,
    # and hence its lower bound must not be negative.

    fit.results <- robustbase::nlrob(
      formula = form,
      data = df,
      method = "M",
      algorithm = "port",
      control = list(warnOnly = TRUE),
      psi = .Mwgt.psi1("huber", cc = 1.345),
      # psi = .Mwgt.psi1("bisquare", cc=4.6),
      start = c(
        Amplitude = 1.0,
        phase = 0,
        offset = offset.start.value
      ),
      lower = lower.limits,
      upper = upper.limits,
      tol = 1e-04,
      trace = FALSE
    )

    return(fit.results)
  }

  # Fit model to data for each site
  set.seed(4334934)
  by_site <- by_site |>
    mutate(
      model = purrr::map(data, site_model),
      `fitted quantity` = quantity.column.name,
      resids = purrr::map2(data, model, add_residuals),
      preds = purrr::map2(data, model, add_predictions),
      fitted_sin_curve = purrr::map2(data, model, iso_add_sinusoidal_prediction)
    )

  return(by_site)
}


iso_extract_residuals <- function(model) {
  residuals <- tidyr::unnest(model, resids) |> 
    dplyr::select(Site, Date, `fitted quantity`, quantity, resid)
  
  return(residuals)
  
}


iso_extract_RMSE_by_site <- function(model) {
  resids.all <- iso_extract_residuals(model)
  RMSE.by.site <- resids.all |> 
    group_by(Site) |> 
    summarise(RMSE = sqrt(sum(resid^2)/n()))
}


iso_extract_MAD_by_site <- function(model) {
  resids.all <- iso_extract_residuals(model)
  MAD.by.site <- resids.all |> 
    group_by(Site) |> 
    summarise(MAD = median(abs(resid)))
}


iso_extract_R_squared_by_site <- function(model, digits = 2) {
  resids.all <- iso_extract_residuals(model)
  R.squared.by.site <- resids.all |>
    group_by(Site) |> 
    mutate(fitted = quantity  - resid, 
           quantity.mean = mean(quantity)) |>
    summarise(R2 = 1 - sum(resid^2)/sum((quantity - quantity.mean)^2)) |>
    mutate(R2 = round(R2, digits))
  
  return(R.squared.by.site)
}


iso_extract_n_samples_by_site <- function(model) {
  resids.all <- iso_extract_residuals(model)
  n.samples.by.site <- resids.all |>
    group_by(Site) |> 
    summarise(n  = n())
  
  return(n.samples.by.site)
}



iso_extract_month_duration <- function(model) {
  # min.date = "2021-01-15"
  # max.date = "2023-08-25"
  # num.months = 31.3
  
  resids.all <- iso_extract_residuals(model)
  Date.table <- resids.all |>
    group_by(Site) |>
    summarise(
      min.date = min(Date),
      max.date = max(Date),
      months.duration.partial = round(interval(min.date, max.date) %/% days(1) /(365/12), 1),
      months.duration =         interval(min.date, max.date) %/% months(1)
    )
  
  return(Date.table)
}


iso_convert_phase_radians_to_days <- function(mydata, phase.radians.column.name) {
  #
  # fitted sinusoidal curve of the form
  #
  # A*sin(x - phase.radians) + offset
  #
  # phase.radians is the phase in radians
  #    - positive value then sin curve shifts to the right
  #    - negative value then sin curve shifts to the left
  #
  # Append columns to mydata for
  #
  #  (1) the phase in days
  #  (2) location of peak in days for the sin curve (can be negative)
  #  (3) location of peak in days for the sin curves: values 0 to 365
  #
  # Usage:
  #
  # iso_convert_phase_radians_to_days(phase.estimate, estimate)
  #
  # where phase.estimate is data frame, for which the column "estimate" 
  # has the phase in radians. 
  #
  
  output  <- mydata |>
    mutate(
      phase.radians = {{ phase.radians.column.name }}, 
      phase.days = round(phase.radians/(2*pi)*365), 
      peak.location.day = round((phase.radians + pi/2)/(2*pi)*365), 
      peak.location.day.positive = dplyr::if_else(peak.location.day < 0, 
                                            peak.location.day + 365, 
                                            peak.location.day)
    ) |>
    dplyr::select(-phase.radians)

  return(output)
}



iso_extract_parameters <- function(model, estimates.only = FALSE) {
  
  summary.table <-  model |>
    mutate(parameters = map(model, broom::tidy)) |>
    unnest(parameters) |> 
    dplyr::select(Site, `fitted quantity`, term, estimate, std.error, statistic, p.value)
  
  if(estimates.only) {
    summary.table <- summary.table |> 
      dplyr::select(Site, term, estimate) |> 
      pivot_wider(names_from =term, values_from = estimate)  
  }
  
  summary.table <- ungroup(summary.table)
  
  return(summary.table)
  
}





iso_add_sinusoidal_prediction <- function(data, model) {
  
  min.date <- min(data$Date)
  max.date <- max(data$Date)
  date.seq <- seq(min.date, max.date, by = "1 day")
  Jdf <- iso_Julian_day_fraction(date.seq)
  
  data.to.predict <- tibble(Date = date.seq, Julian.Day.Fraction = Jdf)
  prediction <- data.to.predict |> 
    modelr::add_predictions(model, var = "preds_fitted_sin_curve")
  
  return(prediction)
  
}



iso_extract_fitted_sin_curve <- function(model) {
  fitted_sin_curve <- model |> 
    dplyr::select(Site, fitted_sin_curve) |> 
    tidyr::unnest(fitted_sin_curve)  

  return(fitted_sin_curve)
  
}



iso_extract_data <- function(model) {
  all.data <- model |> 
    dplyr::select(Site, `fitted quantity`, data) |> 
    tidyr::unnest(data) |> 
    dplyr::select(-Julian.Day.Fraction)
    
  return(all.data)
  
}



iso_plot_raw <- function(data.raw, 
                         site.column.name,
                         date.column.name,
                         quantity.column.name,
                         main.title = NULL,
                         sites = NULL) {
  
  # raw.data: a data frame containing a date column (of Date class) and
  #           other columns:  d18O, d2H 
  #
  # site.column.name: character string for the name of the site id column
  #
  # date.column.name: character string for the name of the date column
  #
  # quantity.column.name: character string for the column to be plotted against the date
  #
  # main.title:  title for plot
  #
  # sites:  vector of characters for selected sites to plot


  data.raw <- data.raw |> 
    dplyr::select(Site = !!site.column.name, 
           Date = !!date.column.name, 
           quantity = !!quantity.column.name)   
  
  # Plot just selected sites
  if(!is.null(sites)) {
    data.raw  <- filter(data.raw,  Site %in% sites)
  }
  
  if(is.null(main.title)) {
    main.title <- paste0("Raw ", quantity.column.name, " by site")
  }
  
  min.year <- min(year(data.raw$Date))
  max.year <- max(year(data.raw$Date))
  min.date <- make_date(year = min.year, month =  1, day =  1)
  max.date <- make_date(year = max.year, month = 12, day = 31)
  seq.years <- seq(from = min.date, to = max.date, by = "years")
  
  ggplot(data.raw, aes(x = Date, y = quantity)) +
    scale_x_date(limits = c(min.date, max.date)) +  
    geom_vline(xintercept = seq.years, linetype = "dashed", colour = "orange") +  
    geom_point(colour = "blue") +
    geom_smooth(colour = "brown") +
    ggtitle(main.title) +
    facet_wrap(vars(Site), ncol = 3, scales = "free_y") 
}



iso_plot_residuals <- function(model, 
                               sites = NULL, 
                               combined = FALSE) {
  residual.values <- iso_extract_residuals(model)
  quantity.text <- unique(residual.values$`fitted quantity`)
  
  if(!is.null(sites)) {
    residual.values <- filter(residual.values, Site %in% sites)
  }
  
  p <- ggplot(residual.values, aes(Date, resid)) +
    geom_point(aes(group = Site), alpha = 1 / 3, colour = "brown") + 
    geom_smooth(se = FALSE, linewidth = 0.5) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    ylab("Residuals (observed - predicted)") +
    ggtitle(paste0(quantity.text, ": residuals")) +
    theme(axis.text.x = element_text(size = 7)) 
  
  if(!combined) p <- p + facet_wrap(vars(Site)) 
  
  p 
  
}



iso_plot_fitted_curve <- function(fitted.models, sites = NULL) {
  
  observed.values  <- iso_extract_data(fitted.models)
  fitted.sin.curve <- iso_extract_fitted_sin_curve(fitted.models)
  quantity.text <- unique(observed.values$`fitted quantity`)
  
  if(!is.null(sites)) {
    observed.values  <- filter(observed.values,  Site %in% sites)
    fitted.sin.curve <- filter(fitted.sin.curve, Site %in% sites)
  }
  
  min.year <- min(year(observed.values$Date))
  max.year <- max(year(observed.values$Date))
  min.date <- make_date(year = min.year, month =  1, day =  1)
  max.date <- make_date(year = max.year, month = 12, day = 31)
  seq.years <- seq(from = min.date, to = max.date, by = "years")  
  
  ggplot(data = observed.values, aes(x = Date, y = quantity)) +
    geom_vline(xintercept = seq.years, linetype = "dashed", colour = "orange") +
    geom_point(colour = "brown", alpha = 1/3) +
    geom_line(data = fitted.sin.curve, 
              aes(Date, preds_fitted_sin_curve),
              colour = "red") +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    ylab(quantity.text) + 
    ggtitle(paste0(quantity.text, ": observed values and fitted curves")) +
    theme_grey(base_size = 8) +
    facet_wrap(vars(Site), scales = "free_y")  
}


iso_plot_fitted_parameters <- function(model, parameter.name, main.title = NULL) {
  
  if(is.null(main.title)) {
    main.title <- paste0("Parameter estimates by site: ", parameter.name)
  }
  
  summary.table <-  model |>
    mutate(parameters = map(model, broom::tidy)) |>
    unnest(parameters) |> 
    dplyr::select(model = Site, `fitted quantity`, term, estimate, std.error, statistic, p.value) |> 
    filter(term == parameter.name)
  
    dotwhisker::secret_weapon(summary.table, 
              var = parameter.name, 
              vline = geom_vline(xintercept = 0, colour = "grey50", linetype = 2)) +
    ggtitle(main.title, subtitle = "(95% confidence interval)")

}



iso_plot_spatial_dynamic <- function(the.data, 
                             site.column.name, 
                             lat.column.name, 
                             lon.column.name,
                             quantity.column.name) {
  
  the.data <- the.data |> 
    dplyr::select(Site = !!site.column.name, 
           lat =  !!lat.column.name, 
           lng =  !!lon.column.name,
           quantity = !!quantity.column.name)
  
  min.quan <- min(the.data$quantity, na.rm = TRUE)
  max.quan <- max(the.data$quantity, na.rm = TRUE)
  
  pal <- leaflet::colorNumeric(c("blue", "green", "red"), c(min.quan, max.quan))
  
  leafMap <- leaflet(data = the.data) |> 
    setView(lat = -41.78, lng = 173.467, zoom = 5) |> 
    addTiles() |>
    addCircleMarkers(color = ~pal(quantity),
                     fillOpacity = 1.0,
                     radius = 5, 
                     # popup = ~htmlEscape(round(quantity,2)), 
                     popup = paste0("Site:", the.data$Site, 
                                    "<br>", 
                                    quantity.column.name, ": ", round(the.data$quantity, 2)),
                     label = ~Site) |> 
    addLegend(pal = pal, 
              values = ~quantity, 
              position = "bottomright",
              className = "legendbox",
              title = quantity.column.name)
  
  leafMap 
}


iso_return_base_nz_map <- function(location = "NZ") {
  
  bbox <- case_when(
    location == "NZ" ~ c(164.81, -47.74, 179.44, -33.79),
    location == "NI" ~ c(171.75, -41.90, 179.44, -33.79),
    location == "SI" ~ c(164.81, -47.74, 174.49, -40.18)
  )
  
  names(bbox) <- c("left", "bottom", "right", "top")  
  
  if(length(bbox) == 1) stop('For the location argument use: "NZ", "NI", "SI"')
  

  nz.base <- ggmap::get_stadiamap(bbox, maptype = "stamen_terrain_background", zoom = 7)

}


iso_plot_site_locations <- function(the.data, 
                                    site.column.name, 
                                    lat.column.name, 
                                    lon.column.name, 
                                    island, 
                                    max.overlaps = 12) {
  
  the.data <- the.data |> 
    dplyr::select(Site = !!site.column.name, 
           lat =  !!lat.column.name, 
           lng =  !!lon.column.name) 
  
 nz.base <- iso_return_base_nz_map(location = island)
 
 ggmap::ggmap(nz.base) +
   geom_point(data = the.data,
              aes(x = lng, y = lat), 
              colour = "red") +
   ggrepel::geom_label_repel(data = the.data, 
              aes(x = lng, y = lat, label = Site), 
              size = 2.5, 
              min.segment.length = 0.5,
              max.overlaps = max.overlaps) +
   xlab("Longitude") +
   ylab("Latitude") +
   ggtitle("Location of sites")

}


iso_plot_spatial_static <- function(the.data, 
                                     site.column.name, 
                                     lat.column.name, 
                                     lon.column.name,
                                     quantity.column.name) {
  
  # 4th October 2023 changed zoom from 7 to 4 due to error:
  # Error in f(init, x[[i]]) : 
  #   number of columns of matrices must match (see arg 2)  
  
  the.data <- the.data |> 
    dplyr::select(Site = !!site.column.name, 
           lat =  !!lat.column.name, 
           lng =  !!lon.column.name,
           quantity = !!quantity.column.name)
  
  bbox <- c(164.81, -47.74, 179.44, -33.79)
  names(bbox) <- c("left", "bottom", "right", "top")
  
  nz.base <- ggmap::get_stadiamap(bbox, maptype = "stamen_terrain_background", zoom = 7) 
 
  
  # Allocate size and colour of circle according to quantity value
  #   - binned for colour
  ggmap::ggmap(nz.base) +
    geom_point(data = the.data, 
               aes(x = lng, y = lat, 
                   size = quantity,
                   colour = quantity), 
    ) +
    scale_color_binned(type = "viridis") +
    labs(colour = quantity.column.name, size = quantity.column.name) +    
    xlab("Longitude") +
    ylab("Latitude") +
    ggtitle(quantity.column.name)
}


iso_fitted_parameters_residuals <- function(model.data, model.fitted) {
  
  predict.params <- predict(model.fitted, newdata = model.data)
  
  predict.params <- as_tibble(predict.params) |> 
    rename(Amplitude.predict = Amplitude, 
           phase.predict = phase, 
           offset.predict = offset) |>   
    mutate(Site = model.data$Site) |> 
    relocate(Site)
  
  model.data <- model.data |> 
    left_join(predict.params)
  
  model.data <- model.data |> 
    mutate(Amplitude.residual = Amplitude - Amplitude.predict,
           Amplitude.residual.percent = Amplitude.residual/Amplitude*100,
           phase.residual = phase - phase.predict, 
           phase.residual.percent = phase.residual/phase*100,
           offset.residual = offset - offset.predict, 
           offset.residual.percent = offset.residual/offset*100)
  
  return.data <- model.data |> 
    dplyr::select(Site, lat, lon, starts_with(c("Amplitude", "phase", "offset")))
  
}




iso_plot_spatial_residual <- function(the.data, 
                                      site.column.name, 
                                      lat.column.name, 
                                      lon.column.name,
                                      quantity.column.name) {
  
  # 4th October 2023 changed zoom from 7 to 4 due to error:
  # Error in f(init, x[[i]]) : 
  #   number of columns of matrices must match (see arg 2)  
  
  the.data <- the.data |> 
    dplyr::select(Site = !!site.column.name, 
           lat =  !!lat.column.name, 
           lng =  !!lon.column.name,
           quantity = !!quantity.column.name)
  
  bbox <- c(164.81, -47.74, 179.44, -33.79)
  names(bbox) <- c("left", "bottom", "right", "top")
  
  nz.base <- ggmap::get_stadiamap(bbox, maptype = "stamen_terrain_background", zoom = 7) 
  
  
  # Allocate size and colour of circle according to quantity value
  #   - binned for colour
  
  the.data <- the.data |> 
    mutate(Sign = ifelse(quantity < 0, "Negative", "Positive"))
  pal <- c("Negative" = "blue", "Positive" = "red")
  
  ggmap::ggmap(nz.base) +
    geom_point(data = the.data,aes(x = lng, y = lat, 
                                   size = abs(quantity), 
                                   colour = Sign)) +
    scale_colour_manual(values = pal, limits = names(pal)) +
    # labs(colour = quantity.column.name, size = quantity.column.name) + 
    labs(size = quantity.column.name) +    
    xlab("Longitude") +
    ylab("Latitude") +
    ggtitle(quantity.column.name)
}



iso_find_catchment_code_for_position <- function(lat, lon, path.to.catchment.files) {
  catchment.code <- NA
  position <- data.frame(long = lon, lat = lat)
  position.sf <- sf::st_as_sf(position, coords = c("long", "lat"),
                              remove = FALSE,
                              crs = 4326)  
  
  all.layers <- sf::st_layers(dsn = path.to.catchment.files)
  
  for (catchment in all.layers$name) {
    the.catchment <- sf::st_read(dsn = path.to.catchment.files, 
                                 layer = catchment, quiet = TRUE)
    within.flag <- position.sf |>
      sf::st_transform(st_crs(the.catchment)) |>
      sf::st_within(the.catchment, sparse = FALSE)
    if(within.flag) {
      catchment.code <- as.numeric(substring(catchment, first = 4))
      break
    }
  }
  
  return(catchment.code)
  
}


iso_find_catchment <- function(locations, 
                               lon.column.name = "lon",                                
                               lat.column.name = "lat", 
                               path.to.catchment.files = "") {
  # locations: data frame contains rows of lat/long values
  # lon.column.name: string with the name of the lon column
  # lat.column.name: string with the name of the lat column

  # path.to.catchment.files: string with path to the catchment shape files
  #
  # Value: number vector of nzreach codes e.g. 2001653
  #
  #
  # locations <- tibble(long = 174.6, lat = -36.3)
  # iso_find_catchment(locations, lat.column.name = "lat", 
  #                               lon.column.name = "long", 
  #                               path.to.catchment.files = "raw data/site_boundaries")
  #
  #
  
  if(path.to.catchment.files == "") { 
    stop("You need to provide an argument for path.catchment.files\n")
  }
  
  the.data <- locations |> 
    dplyr::select(lat =  !!lat.column.name, 
                  lon =  !!lon.column.name)

  # pull out catchments e.g. wsd10010534 
  all.layers <- sf::st_layers(dsn = path.to.catchment.files)
  all.catchments <- all.layers$name

}



iso_find_nzreach_code_for_points <- function(vcsn.lat.lon, 
                                        path.to.catchment.files = "") {
  
  #
  # vcsn.lat.lon:  sf object where each row is a point 
  # path.catchment.file: directory for the catchment shape files
  #
  # Return value: vector of nzreach code for each row
  #
  
  
  if(path.to.catchment.files == "") { 
    stop("You need to provide an argument for path.catchment.files\n")
  }
  
  vcsn.lat.lon$nzreach <- NA
  
   # pull out all catchment labels e.g. wsd10010534 = wsd<nzreach> 
  all.layers <- sf::st_layers(dsn = path.to.catchment.files)
  all.catchments <- all.layers$name

  # Ensure the same CRS as the catchment sf data
  target.crs <- sf::st_crs(all.layers$crs[[1]])  
  vcsn.lat.lon.trans <-  sf::st_transform(vcsn.lat.lon, target.crs)
  
  # Loop over catchment layers and find VCSN points within them  
  for (this.catchment in all.catchments) {
    this.layer <- sf::st_read(dsn = path.to.catchment.files, 
                              layer = this.catchment, 
                              quiet = TRUE)
    # For each row: TRUE within layer, FALSE otherwise
    within.layer <- sf::st_intersects(vcsn.lat.lon.trans, this.layer, sparse = FALSE)
    # "wsd2001653" -> 2001653
    vcsn.lat.lon.trans$nzreach[within.layer] <- as.numeric(substring(this.catchment, first = 4))

  }
  
  # a number e.g. 2001653
  return(vcsn.lat.lon.trans$nzreach)

}


iso_simulate_CI_ratio_truncated_normal <- function(mean1, sd1, mean2, sd2) {
  #
  # Assumed two distributions: 
  #
  # X1 dist normal(mean = mean1, sd = sd1)
  # X2 dist normal(mean = mean2, sd = sd2) 
  #
  # But neither X1 or X2 can be zero or less
  # 
  # Simulate to find 95% confidence interval for X1/X2
  
  num.draws <- 5000
  
  low.bound.cutoff <- mean1 - 1.96*sd1
  low.bound.cutoff <- ifelse(low.bound.cutoff < 0.05, 0.05, low.bound.cutoff)
  numerator   <- rnorm(n = num.draws, mean = mean1, sd = sd1)
  numerator <- ifelse(numerator < low.bound.cutoff,
                      rnorm(n = num.draws, mean = mean1, sd = sd1),
                      numerator)
  
  denominator <- rnorm(n = num.draws, mean = mean2, sd = sd2)

  # 95% confidence interval
  ratio <- numerator/denominator
  ratio <- ifelse(ratio < 0, 0, ratio)
  ratio <- ifelse(ratio > 1, 1, ratio)

  CIs <- quantile(ratio, probs = c(0.025, 0.50, 0.975))
  
  return(list(numerator = numerator, denominator = denominator,
              CIs = CIs, 
              ratio = ratio))
  
}  




iso_low_high_ratio_intervals <- function(mean1, sd1, mean2, sd2) {
  #
  # Assumed two distributions: 
  #
  # X1 dist normal(mean = mean1, sd = sd1)
  # X2 dist normal(mean = mean2, sd = sd2) 
  #
  # But neither X1 or X2 can be zero or less
  # 

  
  
  # lower values for YWF 
  river.low.amp   <- mean1 - 1.96*sd1
  river.low.amp <- ifelse(river.low.amp < 0, 0, river.low.amp)
  river.hgh.amp   <- mean1 + 1.96*sd1
  
  precip.low.amp <- mean2 - 1.96*sd2
  precip.hgh.amp <- mean2 + 1.96*sd2
  
  YWF.low <- river.low.amp/precip.hgh.amp
  YWF.hgh <- river.hgh.amp/precip.low.amp 
  YWF.hgh <- ifelse(YWF.hgh > 1, 1, YWF.hgh)
  
  return(list(YWF.low = YWF.low, YWF.hgh = YWF.hgh))
  
}  



  



