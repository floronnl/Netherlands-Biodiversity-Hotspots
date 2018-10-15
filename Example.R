########################################################
### Example of how to use the HotspotStatus.R script ###
########################################################
#
#
# The HotspotStatus.R script can be used to locate hotspots of biodiversity.
# The script requires unique species observations per grid cell, with optionally a list of species divided in groups.
# Hotspots will be determined per group.
# It is recommended to present the data in a specific way.
# We'll explore a sample dataset provided in the Sample Data folder, it is recommended to order your data in a similar way.
#
# Firstly, we'll have a look at a file containing all grid cell numbers.
#
Gridcells = read.csv("Sample Data\\Gridcells.csv", header = F)$V1
head(Gridcells)
# [1] 13370 13371 14366 14367 14368 14369
#
# Gridcells is a vector of all grid cells we want to study. It is the X and Y coordinates of each cell appended.
# (so X: 40 and Y: 352 becomes grid cell 40352).
#
# Next we'll add a file containing species observations.
#
obs_plants = read.csv("Sample Data\\obs_plants.csv", header = T, sep = ";")
head(obs_plants)
# species  x   y
#        1 14 371
#        1 14 372
#        1 14 377
#        1 14 378
#        1 15 367
#        1 15 370
#
# obs_plants is a data frame of all species observations. The column "species" gives the unique species ID,
# the column "x" is the X coordinate and the column "y" is the Y coordinate of each grid cell.
# Double observations (the same species occurring multiple times in the same grid cell) have already been removed.
#
# Lastly, we'll add a file containing the accepted species and their groups.
# These groups can be functional groups, taxonomical groups etc.
#
plant_species = read.csv("Sample Data\\plant_species.csv", header = T, sep = ";")
head(plant_species)
# accepted_species  group
#                 1 Group1
#                 2 Group2
#                 3 Group1
#                 4 Group2
#                 5 Group1
#                 7 Group2
#
# plant_species is a data frame listing all species we want to take into account for determining the biodiversity hotspots.
# If a species does not occur in this list, but does occur in the observation data frame, it is not used for determining the hotspots.
# The "accepted_species" column gives the unique ID of each species, the column "group" gives the species group of each species.
# Hotspots are determined for each group separately.
#
# To determine which grid cells are hotspots of biodiversity, we first need to load the script.
#
source("HotspotStatus.R")
#
# Now, we can simply execute the "HotspotStatus()" function.
# This will return a table with the results.
# We'll go into each parameter of this function and what it does as well as the results table.
# A description of each parameter can also be found in the HotspotStatus.R script.
#
hotspotResult = HotspotStatus(species = species,               # species is a vector of observed species IDs
                              x = x,                           # x is a vector of X coordinates of the grid cell where the species was observed
                              y = y,                           # y is a vector of Y coordinates of the grid cell where the species was observed
                              observationTable = obs_plants,   # Optional. Data frame containing the above three vectors. The vectors can also be provided separately.
                              speciesList = accepted_species,  # Optional. Vector of species IDs to be taken into account for determining the hotspots. If not provided all observed species will be considered.                
                              speciesGroups = group,           # Optional. Vector of groups the species are part of. Hotspots will be determined for each group. If not provided, all species will be considered to be in one group.
                              speciesTable = plant_species,    # Optional. Data frame containing the above two vectors. The vectors can also be provided separately.
                              gridcellList = Gridcells,        # Vector of all grid cells. Grid cells need to be the X and Y coordinate appended.
                              numberOfHotspots = 2,            # The number of hotspots per group to be determined. If there are ties, the actual number of hotspots might be slightly higher.
                              subTopMode = "amount",           # If and how the subtop should be calculated ("none", "amount", "sd", "2sd"). This subtop is used to extend the hotspot to be larger than a single grid cell.
                              subTop = 15,                     # Number of subtop cells if subTopMode is "amount".
                              returnTF = T,                    # Whether the results should print TRUE/FALSE (T) or 1/0 (F).
                              printProgress = F                # Whether progress should be printed. Will slightly slow down the script, but might be nice with large datasets.
                              )
# The script will now do the following things:
# 1. Determine the amount of species in each grid cell.
# 2. Determine for each species group separately whether each grid cell belongs in its subtop (in this case, 15 subtop grid cells will be determined per group).
# 3. For each species group separately, determine the grid cell with the single highest biodiversity.
# 4. Find subtop grid cells adjacent to this cell, adding these to the hotspot.
# 5. Find the next grid cell with the highest biodiversity that is not yet part of a hotspot.
# 6. Repeat steps 4 and 5 until the desired number of hotspots has been found (in this case, 2 per species group).
#
# The results table has been saved in "hotspotResult", let's have a look at it.
#
head(hotspotResult)
#  Gridcell  x   y Group1 Hotspot_Group1 subtop_Group1 hotspot_ID_Group1 Group2 Hotspot_Group2 subtop_Group2 hotspot_ID_Group2
#     13370 13 370      0          FALSE         FALSE                 0      1          FALSE         FALSE                 0
#     13371 13 371      5          FALSE         FALSE                 0      4          FALSE         FALSE                 0
#     14366 14 366      0          FALSE         FALSE                 0      0          FALSE         FALSE                 0
#     14367 14 367      1          FALSE         FALSE                 0      1          FALSE         FALSE                 0
#     14368 14 368      0          FALSE         FALSE                 0      0          FALSE         FALSE                 0
#     14369 14 369      1          FALSE         FALSE                 0      2          FALSE         FALSE                 0
#
# We'll now go over each column in the results.
#
# The "Gridcell" column contains each grid cell that was investigated.
# The "x" column contains the X coordinate of each grid cell.
# The "y" column contains the Y coordinate of each grid cell.
# The next couple of columns will be repeated for each species group.
# The column with the name of the species group ("Group1" in this case) contains the number of species of that group found in each grid cell.
# The column "Hotspot_" followed by the group name ("Hotspot_Group1" in this case) contains whether each grid cell is the centre of a hotspot for that species group.
# The column "subtop_" followed by the group name ("subtop_Group1" in this case) contains whether each grid cell is within the subtop for that species group. Note: Subtop cells are not necessarily in a hotspot!
# The column "hotspot_ID_" followed by the group name ("hotspot_ID_Group1" in this case) contains the hotspot ID for each grid cell that group.
# If the hotspot ID is 0, the grid cell is not part of a hotspot for that species group.
#
# We now want to only view the results for the grid cells that are part of a hotspot for species group "Group1":
#
hotspotResult[hotspotResult$hotspot_ID_Group1 != 0 ,]
# Gridcell  x   y Group1 Hotspot_Group1 subtop_Group1 hotspot_ID_Group1 Group2 Hotspot_Group2 subtop_Group2 hotspot_ID_Group2
#     22380 22 380     11           TRUE          TRUE                 1      8          FALSE         FALSE                 0
#     25388 25 388      8          FALSE          TRUE                 2      7          FALSE         FALSE                 0
#     25389 25 389     10           TRUE          TRUE                 2     10          FALSE          TRUE                 0
#     25399 25 399     10           TRUE          TRUE                 3     10          FALSE          TRUE                 0
#     26400 26 400      8          FALSE          TRUE                 3     10          FALSE          TRUE                 0
#     27380 27 380     10           TRUE          TRUE                 4     10          FALSE          TRUE                 0
#     27400 27 400      9          FALSE          TRUE                 3     10          FALSE          TRUE                 0
#     31391 31 391      8          FALSE          TRUE                 5     11           TRUE          TRUE                 3
#     32391 32 391     10           TRUE          TRUE                 5     10          FALSE          TRUE                 3
#
# As you can see, 5 hotspots have been found for "Group1", although we said we only wanted 2.
# This is because there were 4 grid cells with 10 species within "Group1".
# Hotspots with more than 1 cell with the same ID are larger than 1 grid cell.
# It is recommended to visualize the hotspots on a map using GIS software. This can be done by giving each grid cell a colour dependent on the hotspot ID.
#
##