#install.packages("ggplot2")
#install.packages("ggrepel")
#install.packages("readstata13")

library(ggplot2)
library(ggrepel)

# setwd("D:\\DATA")
# fnm = 'VOLCANOPLOT_DATASET.dta'

fnmCSV = '/zdsk/Manuscripts/Baydoun/UK-BioBank/2024-02-13--13--ProteomicCVD/Volcano/2024-10-18/LE8_PROTEOME_DMRIVOLCANOPLOTDATA.csv'
datObj = read.csv(fnmCSV, as.is=T)
nrow(datObj)
summary(datObj)
head(datObj)

pltDat = within(datObj, {
	lblOutcome	= sub('^ z', '', Outcome, perl=T)
# 	p		= pmax(p, p, 10^-600)
	})

summary(pltDat)
head(pltDat)

options(ggrepel.max.overlaps = Inf)	# avoids error message about colliding labels

Volcano = function(volDat) {
	ggObj =	ggplot(data=volDat, aes(x=estimate, y=mlog10p, label=lblOutcome)) +
		geom_point(shape=20, size=1) +
                geom_point(shape=20, size=3, col='red',    data=subset(volDat, selected_final==1)) +
		geom_point(shape=20, size=3, col='orange', data=subset(volDat, selected_final==1 & signif==1 & estimate>  0.20)) +
                geom_point(shape=20, size=3, col='blue',   data=subset(volDat, selected_final==1 & signif==1 & estimate< -0.20)) +
		xlim(-.45, .45) + ylim(0, 250) +
		geom_text_repel(size=3, col='red', data=subset(volDat,  selected_final==1 & estimate>  0.20), aes(label=lblOutcome)) +
                geom_text_repel(size=3, col='blue', data=subset(volDat, selected_final==1 & estimate< -0.20), aes(label=lblOutcome)) +
                labs(x='Effect size') +
		theme_minimal() +
		geom_vline(xintercept=-.20, col="darkgray", linetype=2) +
		geom_vline(xintercept=.20,  col="darkgray", linetype=2) +
		geom_hline(yintercept=-log10(0.05/1463), col="darkgray", linetype=2)

	ggObj
	}

jnk = grid::current.viewport()	#prevents UseMethod("depth") error

(pltObj = Volcano(pltDat))

dirName = '/zdsk/Manuscripts/Baydoun/UK-BioBank/2024-02-13--13--ProteomicCVD/Volcano/2024-10-18/'
imgName = '2024-10-18-UKB13-Volcano'
pdfFile = paste0(dirName, imgName, '.pdf')
pngFile = paste0(dirName, imgName, '.png')

pdf(pdfFile)
pltObj
jnk = dev.off()

png(pngFile)
pltObj
jnk = dev.off()
