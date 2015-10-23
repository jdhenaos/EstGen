library("GWASTools")

SNPs<-read.table("SNPsTaller4.txt")
genotype<-as.matrix(SNPs[2:16,2:5])
SNPID <- 1:15
chrom<-as.integer(c(10,3,6,19,6,6,23,18,11,23,1,13,6,12,1))
pos<-as.integer(c(6110829,28071444,159465977,6668972,135739355,137452908,34892503,28789725,3259809,13615118,182549019,50185204,29705659,6440009,117104215))
IDscan<-1:4
InputData<-MatrixGenotypeReader(genotype=genotype,snpID=SNPID,chromosome=chrom,position=pos,scanID=IDscan)


sex <- c('F','F','F','F')

sex

scanID <- as.numeric(IDscan)
class(scanID)
frame <- data.frame(scanID,sex)

scan <- ScanAnnotationDataFrame(frame)

getScanID(scan)
getSex(scan)


InputData2<-GenotypeData(InputData, snpAnnot = NULL, scanAnnot = scan)


fisher <- batchFisherTest(InputData2, batchVar = "scanID",
                          chrom.include=23, sex.include = "F",return.by.snp=TRUE)
fisher
