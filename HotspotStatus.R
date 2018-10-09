#####################################
#### Created by Dion van der Hak ####
#####################################
#
# Version: 0.5.0
#
### Description: ###
#
# This function takes the observations as input, as well as the species list and a list of all gridcells.
# There must not be any double observations (same species in same gridcell).
# The species list must only contain species that are counted, if a species was observed but does not occur in the
# species list, it is ignored.
# If no species list is given, all species are counted.
# The gridcell list must contain all relevant gridcells, also ones without observations.
# The gridcell names must be unique and the x and y coordinate appended (so x: 40 and y: 352 becomes gridcell: 40352).
#
# This function returns a dataframe with for each gridcell whether or not they meet the Hotspot requirement.
#

HotspotStatus = function(species,                             # Vector of observed species
                         x,                                   # Vecor of x coordinates of the observations
                         y,                                   # Vecor of y coordinates of the observations
                         observationTable = NULL,             # Not required, dataframe containing above 3 vectors
                         speciesList = NULL,                  # Vector of species to be counted, if not present all species will be counted
                         speciesGroups = c("Biodiversity"),   # Not required, vector of species groups the species in speciesList belong in
                         speciesTable = NULL,                 # Not required, dataframe containing above 2 vectors
                         gridcellList,                        # Vector of all gridcells
                         numberOfHotspots = 10,               # Number of hotspots to be found
                         subTopMode = "none",                 # How the subtop should be calculated ("none", "amount", "sd", "2sd")
                         subTop = 250,                        # Number of subtop gridcells (if subTopMode is amount)
                         returnTF = T,                        # Whether hotspot status should be returned as true/false or 1/0
                         printProgress = F                    # Whether calculation progress should be printed
                         ) {                      

  
  ## Init ##
  # Returns error if one of the entered variables has the wrong datatype
  stopifnot(is.numeric(numberOfHotspots))
  stopifnot(is.numeric(subTop))
  stopifnot(is.logical(returnTF))
  stopifnot(is.logical(printProgress))
  
  # Possibilities for subTopMode, disables subTopMode if not allowed
  subTopModeAllowed = c("none", "amount", "sd", "2sd")
  if(!subTopMode %in% subTopModeAllowed) {
    subTopMode = "none"
  }

  # Load the required libraries
  require(reshape2)
  require(plyr)
  
  # Fix variables to allow for the user to point to the dataframe
  species = eval(substitute(species), observationTable, parent.frame())
  x = eval(substitute(x), observationTable, parent.frame())
  y = eval(substitute(y), observationTable, parent.frame())
  speciesList = eval(substitute(speciesList), speciesTable, parent.frame())
  speciesGroups = as.character(eval(substitute(speciesGroups), speciesTable, parent.frame()))
  
  # Create gridcell variable
  Gridcell = as.numeric(paste0(x, y))
  
  # Generate species list if absent
  if(is.null(speciesList)) {
    speciesList = unique(species)
  }
  
  
  ## Functions ##
  # Get x coordinate from gridcell number
  getx = function(gridcell) {
    return(as.numeric(substr(as.character(gridcell), 1, nchar(gridcell)-3)))
  }
  
  # Get y coordinate from gridcell number
  gety = function(gridcell) {
    return(as.numeric(substr(as.character(gridcell),nchar(gridcell)-2,nchar(gridcell))))
  }
  
  # Get the required richness for the Hotspot requirement
  getHotspotRequirement = function(richness, numberOfHotspots) {
    requiredRichness = min(head(sort(richness[,1], decreasing = T), numberOfHotspots))
    return(requiredRichness)
  }
  
  # Test whether each grid cell meets the Hotspot requirement
  testHotspotStatus = function(requiredRichness, richnessVector) {
    Hotspot_status = richnessVector >= requiredRichness
    return(Hotspot_status)
  }
  
  # Checks the surroundings of a given gridcell for subtop gridcells
  # gridcell is a list of all gridcells, of which x and y are coordinates
  # Hotspot is the gridcell number of the hotspot center: the cell whose surroundings need to be checked
  # subtop is whether gridcell is a subtop or not (T/F)
  checkSurroundings = function(gridcell, x, y, Hotspot, subtop, data = NULL) {
    # Allow for data substitution
    gridcell = eval(substitute(gridcell), data, parent.frame())
    x = eval(substitute(x), data, parent.frame())
    y = eval(substitute(y), data, parent.frame())
    subtop = eval(substitute(subtop), data, parent.frame())
    # Bind
    data = data.frame(gridcell, x, y, subtop)
    # Check
    data$newCell = "" # New cells are cells whose surroundings haven't been checked yet
    data$newCell[data$gridcell == Hotspot] = T # All Hotspot cells are new
    data$hotspot = F # Create hotspot column
    data$hotspot[data$gridcell == Hotspot] = T # Hotspot cells are hotspots

    while(length(data$newCell[data$newCell == T]) != 0) { # Keep running until no more new hotspot gridcells have been found
      for(i in data$gridcell[data$newCell == T]) { # For each hotspot cell whose surroundings haven't been checked yet
        data$newCell[data$gridcell == i] = F # Cell is no longer new since it's surroundings are going to be checked
        x = data$x[data$gridcell == i]
        y = data$y[data$gridcell == i]
        gridcell.surround = c(paste0(x-1, y-1), paste0(x, y-1), paste0(x+1, y-1),
                              paste0(x-1, y), paste0(x+1, y),
                              paste0(x-1, y+1), paste0(x, y+1), paste0(x+1, y+1)) # Define coordinates of surrounding cells
        for(j in gridcell.surround) { # For each surrounding cell
          if(!j %in% data$gridcell) {
            next # Stop if cell doesn't exist
          }
          if(data$subtop[data$gridcell == j] == T) {
            data$hotspot[data$gridcell == j] = T # If the surrounding cell is in the subtop, set it as a hotspot cell
            if(data$newCell[data$gridcell == j] != F) {
              data$newCell[data$gridcell == j] = T # Also: if the hotspot is new, set it as such
            }
          }
        }
      }
    }
    return(data$gridcell[data$hotspot == T]) # Returns vector of gridcells that are part of a hotspot
  }
  

  ## Main ##
  # Bind observations together into table
  observationTable = data.frame(species, Gridcell, x, y)
  
  # Bind species together into table
  speciesTable = data.frame(speciesList, speciesGroups)
  
  # Transform the gridcell list into a table
  gridcellList = cbind(gridcellList)
  colnames(gridcellList) = c("Gridcell")
  
  # Remove observations of species that do not occur in the species list
  observationTable = observationTable[observationTable$species %in% speciesTable$speciesList,]
  
  # Convert the observations table to richness table
  groups = paste0(unique(speciesTable$speciesGroups)) # Makes variable for the different groups
  richnessTable = data.frame(gridcellList) # Prepares the richness table
  richnessTable$x = getx(richnessTable$Gridcell) # Adds x coordinate
  richnessTable$y = gety(richnessTable$Gridcell) # Adds y coordinate

  
  for(i in groups) { 
    if(printProgress == T) { # Prints progress
      print(paste0("Calculating hotspots for ", i))
    }
    
    # Counts the number of species for each group
    temp = as.data.frame( 
      table(
        observationTable$Gridcell[
          observationTable$species %in% speciesTable$speciesList[
            speciesTable$speciesGroups == i]]))
    colnames(temp) = c("Gridcell", i)
    temp = merge(gridcellList, temp, by = "Gridcell", all = T) # Merges the richness table with the gridcell list
    temp[is.na(temp)] = 0 # Removes NAs
    richnessTable = merge(richnessTable, temp, by = "Gridcell", all = T) # Adds current group to the richness table
    richnessTable[is.na(richnessTable)] = 0 # Removes NAs
    
    # This part is only to allow for subtop calculation of sd and 2sd, might be replaced later
    # Get Hotspot requirement (gets the minimum richness required to meet the Hotspot requirement given numberOfHotspots).
    requiredRichness = max(1, getHotspotRequirement(richnessTable[i], numberOfHotspots))
    # Test Hotspot status
    richnessTable[paste0("Hotspot_", i)] = as.vector(testHotspotStatus(requiredRichness, richnessTable[i]))
    
    # Get subtop requirement if subTopMode is amount (gets the minimum richness required to be subtop given subTop)
    if(subTopMode == "amount") {
      requiredSubTopRichness = max(1, getHotspotRequirement(richnessTable[i], subTop))
    }
    # Get subtop requirement if subTopMode is sd (this is the lowest Hotspot richness minus the sd)
    if(subTopMode == "sd") {
      requiredSubTopRichness = max(1, min(richnessTable[i][richnessTable[paste0("Hotspot_", i)] == T]) - sd(richnessTable[i][richnessTable[paste0("Hotspot_", i)] == T]))
    }
    # Get subtop requirement if subTopMode is 2sd (this is the lowest Hotspot richness minus 2 times sd)
    if(subTopMode == "2sd") {
      requiredSubTopRichness = max(1, min(richnessTable[i][richnessTable[paste0("Hotspot_", i)] == T]) - 2*sd(richnessTable[i][richnessTable[paste0("Hotspot_", i)] == T]))
    }
    # Test subtop status
    if(subTopMode != "none") {
      richnessTable[paste0("subtop_", i)] = as.vector(testHotspotStatus(requiredSubTopRichness, richnessTable[i]))
    }
      
    # Create hotspots
    hotspotList = c() # List of hotspot IDs that have been added
    richnessTable[paste0("Hotspot_", i)] = F # No Hotspots have been defined yet
    richnessTable[paste0("hotspot_ID_", i)] = 0 # All gridcells have hotspot ID 0
    tempNumberOfHotspots = numberOfHotspots
    while(length(hotspotList) < tempNumberOfHotspots) { # New hotspots will be created until the correct amount has been reached
      m = richnessTable$Gridcell[richnessTable[i] == max(richnessTable[i][richnessTable[paste0("hotspot_ID_", i)] == 0]) & richnessTable[paste0("hotspot_ID_", i)] == 0][1] # Gridcell number of cell with the highest richness that is not a hotspot yet
      n = length(hotspotList) + 1 # ID for current hotspot
      richnessTable[paste0("Hotspot_", i)][richnessTable$Gridcell == m,] = T # Sets current Hotspot cell as Hotspot cell
      richnessTable[paste0("hotspot_ID_", i)][richnessTable$Gridcell == m,] = n # Sets hotspot ID for this Hotspot cell
      if(subTopMode != "none") {
        temp = checkSurroundings(Gridcell, x, y, m, get(paste0("subtop_", i)), data = richnessTable) # Check surroundings of each Hotspot cell for subtop cells, returning gridcell numbers of all valid cells
        richnessTable[paste0("hotspot_ID_", i)][richnessTable$Gridcell %in% temp,] = n # Sets hotspot ID for accepted hotspot cells
      }
      hotspotList = append(hotspotList, n) # Makes a list of hotspot IDs
      if(length(hotspotList) >= numberOfHotspots & length(richnessTable$Gridcell[richnessTable[i] == richnessTable[i][richnessTable$Gridcell == m,] & richnessTable[paste0("hotspot_ID_", i)] == 0]) > 0) { # Includes tied last places into the hotspot list
        tempNumberOfHotspots = tempNumberOfHotspots + 1
      }
    }
  }
  
  # Convert all TRUE / FALSE to 1 / 0 if returnTF is false
  if(returnTF == F) {
    richnessTable[richnessTable == F] = 0
    richnessTable[richnessTable == T] = 1
  }
  
  if(printProgress == T) { # Prints done
    print("Done!")
  }
  

  # Return results
  return(richnessTable)
}


### Changelog: ###
#
# 0.5.0: Cleaned up for GitHub publication
# 0.4.1: Fixed crash when subTopMode is "none"
# 0.4.0: Ties will now be included in the hotspot list
# 0.3.0: Added subTopMode "sd" and "2sd"
# 0.2.0: Optimized neighbour search
# 0.1.0: Created function
#
