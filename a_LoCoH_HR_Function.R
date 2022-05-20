library(tidyverse)
library(sf)
library(parallel)

a_LoCoH_HR <- function(data, id = "ID", date = "DATE", 
                       min_days = 30, min_fixes = 30, 
                       iso_levels = c(0.95, 1), fill_holes = TRUE, 
                       min_neighbors = 2, n_cores = detectCores()- 2) {
  
  # Summarise the number of days an individual was tracked and 
  # the number of fixes in the dataset
  data_summary = data %>%
    st_drop_geometry() %>%
    summarise(Days_Tracked = difftime(max(.[[date]], na.rm = TRUE), 
                                      min(.[[date]], na.rm = TRUE), 
                                      units = "days"),
              n_Fixes = n())
  
  # If the number of days or fixes is less than the user defined 
  # threshold (defaults to 30 days/fixes), do not estimate HR and 
  # return number of days/fixes with message
  if(data_summary$Days_Tracked < min_days | 
     data_summary$n_Fixes < min_fixes){
    
    capture.output(paste("Insufficient data to estimate home range:", 
                         unique(data[[id]]),
                         "n Days =",
                         data_summary$Days_Tracked, 
                         "n Fixes = ", data_summary$n_Fixes, sep = " "))
    
    } else{
    
    # Print start time for NN identification to console
    start_time = Sys.time()
    print(paste("Identifying Nearest Neighbors for", 
                nrow(data), "Points:", 
                unique(data[[id]]), Sys.time()))
    
    
    # Create a dataframe with the pairwise distance between 
    # each point in the telemetry data
    pairwise_distances = st_distance(data) %>%
      as_tibble(rownames = "Row_Number", 
                .name_repair = ~paste0("V", .x)) %>%
      mutate(across(everything(), 
                    as.numeric))
    
    # Gets the maximum distance between any two points in 
    # the telemetry data. Will serve as the "a" parameter
    max_distance = pairwise_distances %>%
      select(-Row_Number) %>% max()
    
    # Set up parallelization:
    # If running on linux/mac, will utilize FORK clusters which are 
    # faster to initialize and more memory efficient
    if(Sys.info()[["sysname"]] %in% c("Linux", "Darwin")){
      
      cl <- makeCluster(n_cores, type = "FORK")
      
    } 
    # If running on windows, will utilize SOCKET clusters,
    # which require explicit exporting of any data/packages 
    # used in the function being run
    else{
      cl <- makeCluster(n_cores)
      clusterExport(cl, c("data", 
                          "pairwise_distances", 
                          "min_neighbors",
                          "max_distance"), 
                    envir = environment())
      clusterEvalQ(cl, {
        library(tidyverse)
        library(sf)
      })
    }
    
    # The function below will calculate the convex hull for each 
    # focal point (each row in the dataframe)
    hulls = parLapply(cl,
                      1:nrow(data), 
                      function(row){
                        
                        # From the dataframe with the pairwise distances 
                        # for all points, select the row number and 
                        # column representing the distance to the focal point. 
                        nearest_neighbors = pairwise_distances %>%
                          select(Row_Number, Distance = paste0("V", row)) %>%
                          # Arrange the dataframe by distance
                          arrange(Distance) %>%
                          # Calculate the cumulative distance and filter 
                          # so that only rows where the cumulative distance
                          # is less than "a" (the maximum distance between any 
                          # two points in the dataset) are returned OR 
                          # to meet the specified minimum number of nearest neighbors
                          mutate(Cumulative_Distance = cumsum(Distance)) %>%
                          filter(Cumulative_Distance <= max_distance | 
                                   row_number() <= 1 + min_neighbors) %>%
                          # Get the row number representing each point identified 
                          # as a nearest neighbor
                          pull(Row_Number)
                        
                        # Using those points, create the convex hull 
                        # for the focal point and calculate the number
                        # of points contained within the hull and the 
                        # area of the hull
                        hull = data %>% 
                          filter(row_number() %in% nearest_neighbors) %>%
                          st_union() %>%
                          st_convex_hull() %>%
                          st_as_sf() %>%
                          mutate(n_points = length(nearest_neighbors),
                                 area = st_area(.))
                        
                      }) %>%
      # Bind hulls into a single dataframe and arrange by the 
      # number of points contained within the hull. If there is 
      # a tie in the number of points contained within hulls, 
      # they are arranged by area, from smallest to largest
      bind_rows() %>%
      arrange(desc(n_points), area)
    
    stopCluster(cl)
    
    print(paste("Hulls Created:", unique(data[[id]]), Sys.time()))
    
    # This function will repeat until the number of points contained 
    # within the union of the hulls is greater than the user defined
    # iso_level. For example, if iso_level = 0.95, it will break the loop 
    # when a percent of greater than 95% of points contained within the hull 
    # is reached, and return the largest union of hulls with < 95% of points 
    # contained
    Isopleths = lapply(iso_levels, 
                       function(iso_level){
                         
                         # If the iso_level is 1, return the 100% isopleth, 
                         # otherwise get desired isopleth below
                         if(iso_level == 1){
                           
                           Isopleth = hulls %>% 
                             st_union() %>%
                             st_as_sf() %>%
                             mutate(n_contained = lengths(st_intersects(., data)),
                                    percent = n_contained/nrow(data), 
                                    iso_level = iso_level)
                           
                         }else{
                           
                           # Rather than starting the repeat loop with the hull containing the
                           # the greatest number of points, this creates a sequence ranging from 
                           # 1 to the number of points in the dataset (also the number of hulls). 
                           # Using this sequence, calculate the percent of points contained within 
                           # each hull, then select an initial hull index value (the hull with the 
                           # largest number of points contained that is less than the use defined iso
                           # level).
                           initial_hull_numbers = unique(round(seq(1, nrow(data), length.out = 100)))
                           
                           hull_index = lapply(initial_hull_numbers,
                                               function(hull_index){ 
                                                 
                                                 hulls %>% 
                                                   head(hull_index) %>%
                                                   st_union() %>%
                                                   st_as_sf() %>%
                                                   mutate(n_contained = lengths(st_intersects(., data)),
                                                          percent = n_contained/nrow(data), 
                                                          iso_level = iso_level) 
                                               }) %>%
                             bind_rows(.id = "Index") %>%
                             filter(percent < iso_level) %>%
                             filter(percent == max(percent)) %>%
                             filter(Index == max(Index)) %>%
                             pull(Index)
                           
                           hull_index = initial_hull_numbers[as.numeric(hull_index)]
                           
                           # From the initial hull index identified above, begin the repeat loop
                           
                           repeat{
                             
                             # Create isopleth and test if greater than or
                             # equal to the defined iso_level
                             test_isopleth = hulls %>%
                               head(hull_index) %>%
                               st_union() %>%
                               st_as_sf() %>%
                               mutate(n_contained = lengths(st_intersects(., data)),
                                      percent = n_contained/nrow(data),
                                      iso_level = iso_level)
                             
                             # Break the loop as described above
                             if(test_isopleth$percent >= iso_level){
                               break
                             }
                             
                             # If test_isopleth percent is not greater than iso_level
                             # set as new Isopleth and add 1 to hull_index to restart
                             # the loop
                             hull_index = hull_index + 1
                             
                             Isopleth = test_isopleth
                             
                           }
                           
                         }
                         
                         print(paste(iso_level, "Isopleth Identified:", unique(data[[id]]), Sys.time()))
                         return(Isopleth)
                         
                       })
    
    # If the fill hole argument is TRUE,  
    # fill_holes function with threshold and 
    # recalculate HR area after filling holes 
    if(fill_holes == TRUE){
      
      Isopleths = Isopleths %>%
        map_df(., ~.x %>%
                 mutate(ID = unique(data[[id]])) %>%
                 sfheaders::sf_remove_holes() %>%
                 mutate(Area = st_area(.)) %>%
                 rename(geometry = x) %>%
                 st_transform(4326))
      
    } else{
      
      Isopleths = Isopleths %>%
        map_df(., ~.x %>%
                 mutate(ID = unique(data[[id]]),
                        Area = st_area(.)) %>%
                 rename(geometry = x) %>%
                 st_transform(4326))
      
    }
    
    gc()
    time_difference = difftime(Sys.time(), start_time)
    print(time_difference)
    Isopleths
  }
}
