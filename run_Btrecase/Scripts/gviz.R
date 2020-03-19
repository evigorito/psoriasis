library(rtracklayer)
library(GenomicRanges)
library(ggplot2)

## ebg <- readRDS("/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_ebg.rds")

mySession = browserSession("UCSC")
genome(mySession) <- "hg19"
track.name <-  "gtexGene"

t <- tableNames(ucscTableQuery(mySession, track=track.name))
getTable(ucscTableQuery (mySession, track=track.name, range=e2f3.grange, table=t[1]))

track(session, "targets") <- targetTrack
subTargetTrack <- targetTrack[1] # get first feature
##view <- browserView(session, subTargetTrack * -10, pack = "targetScanS")
loaded_tracks <- trackNames(session)
subTargetTrack <- track(session, "targetScanS")



##from <- 65921878
##to <- 65980988
## knownGenes <- UcscTrack(genome = "mm9", chromosome = "chrX",    track = "xenoRefGene", from = from, to = to,     trackType = "GeneRegionTrack", rstarts = "exonStarts",     rends = "exonEnds", gene = "name", symbol = "name2",     transcript = "name", strand = "strand", fill = "#8282d2",     stacking = "dense", name = "Other RefSeq")

## plotTracks(knownGenes, from=from, to=to, collapseTranscripts = "longest", transcriptAnnotation="symbol")
## z <- ranges(knownGenes)

## pgIslands <- UcscTrack(genome = "mm9", chromosome = "chrX",     track = "cpgIslandExt", from = from, to = to,     trackType = "AnnotationTrack", start = "chromStart",     end = "chromEnd", id = "name", shape = "box",     fill = "#006400", name = "CpG Islands")

######################################################################################################
############# Define GeomGene #######################################################################


## GeomGene = ggproto("GeomGene", Geom,
##                    required_aes = 'x',
##                    default_aes = aes(shape = 19, colour = 'black',
##                                      fill = 'green4', size = 3,
##                                      linetype = 1, alpha = 1,
##                                      fontsize = 1),
##                    draw_key = draw_key_point,

##                    draw_group = function(data, panel_params, coord){

##                        common <- list(
##                            colour = data$colour,
##                            size = data$size,
##                            linetype = data$linetype,
##                            fill = alpha(data$fill, data$alpha),
##                            group = data$group
##                        )

                       

