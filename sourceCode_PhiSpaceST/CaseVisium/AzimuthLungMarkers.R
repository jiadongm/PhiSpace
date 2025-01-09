AzimuthMarkers <- read.csv(
  paste0(dat_dir, "data/LungRef/AzimuthMarkers.csv"),
  header = F
)

dat <- AzimuthMarkers |>
  as.data.frame()

outList <- vector("list", nrow(dat)) |>
  `names<-`(dat[,1])
for(i in 1:nrow(dat)){

  outList[[i]] <- stringr::str_split_1(dat[i,2], ", ")
}

outList <- outList[order(names(outList))]

names(outList) <- gsub("φ", "ph", names(outList))

names(outList)[37] <- "Monocyte-derived Mph"
names(outList)[43] <- "Non-classical monocytes"
names(outList)[57] <- "AT0"

outList <- outList[order(names(outList))]
saveRDS(outList, paste0(dat_dir, "data/LungRef/AzimuthLungMarkers.rds"))

