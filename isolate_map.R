
#### isolate map

## analyses for McDaniel et al done with R version 3.5.3
## ggmap version 2.6.1

install.packages("ggmap")
library("ggmap")

## County Lat Long
## Durham Co  36.02861111  -78.92444444
## Dinwiddie Co 37.10166667 -77.54333333
## Prince Georges Co  38.82944444 -76.87444444
## Adams Co 39.90638889 -77.52972222
## Westchester Co 41.27472222 -73.71055556
## Rensselaer Co  42.72333333 -73.27777778
## Essex Co 44.25083333 -73.81888889


locations.lat <- c(36.02861111,37.10166667,38.82944444,39.90638889,41.27472222,42.72333333,44.25083333)
locations.lon <- c(-78.92444444,-77.54333333,-76.87444444,-77.52972222,-73.71055556,-73.27777778,-73.81888889)

world <- map_data("world")
lakes <- map_data("lakes")
states <- map_data("states")

locations_text <- c("D", "N", "S", "A", "W", "R", "E")

png("figures/isolate_map_2019.png", width = 8, height = 6, units = 'in', res = 300)

mp <- NULL
mp <- ggplot() + geom_polygon(data = world,aes(x=long, y=lat, group = group), 
                              fill="gray65") + coord_fixed(1.3) + borders("state")
mp <- mp + geom_point(aes(x=locations.lon, y=locations.lat),color="black", size=4)
mp <- mp + geom_text(aes(x=locations.lon-1, y=locations.lat, label = locations_text),
                     size=7)
mp <- mp + geom_polygon(data = lakes,aes(x=long, y=lat, group = group), 
               fill="white") + coord_fixed(1.3)
mp + coord_fixed(xlim = c(-82.5, -67.5),ylim = c(34,46), ratio = 1.3) +
  theme(panel.background=element_rect(fill="white"), plot.title=element_text(size=30, hjust=0.5),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=20),
        axis.text=element_text(size=15))


dev.off()
    

