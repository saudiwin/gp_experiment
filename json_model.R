#Load JSON data and model user locations
#Robert Kubinec v0.1

require(jsonlite)

beacon_data <- fromJSON(file('dataset/beaconData.json'))
