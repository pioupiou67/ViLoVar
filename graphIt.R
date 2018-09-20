#!/usr/bin/Rscript --slave

argv <- commandArgs(TRUE)
# on recupere les genes en arguments
strGene <- as.character(argv[1])
# et on les converti en list
listGene <- as.list(strsplit(strGene, ",")[[1]])

graphIt=function(listGene){
  ####### varbar #######
  library(ggplot2)
  #['ULK4','TRPM6','TMEM39A']
  #print(listGene)
  for (g in listGene) {
    print(g)
    ficName <- paste("results/", g, "_R.tsv", sep="")
    t <- read.csv(ficName, header=TRUE, sep = "\t", dec = ",")
    summary(t)
    
    ###couleur/size = effect/freq###
    ####colour by effect
    p2 <- ggplot(data=t, aes(x=position, y=indiv, colour=effect, shape=hasRsId, size=occurVar)) + geom_point()
    #p2 <- ggplot(data=t, aes(x=position, y=indiv, colour=effect, shape=zyg, size=occurVar)) + geom_point()
    ####colour by effect
    p2 <- p2 + scale_colour_manual(values = c("HIGH"= "red","HIGHsplice"="orange","MODERATE"="blue","MODERATEsplice"="lightblue","LOW"="green","LOWsplice"="yellow"))
    
    #### Add colours to the name selon nb de variant dans le gene (colonne 8)
    d=c()
    h=length(t[,1])
    
    if (t[1,7]>1){
      d=c(d,"red")
    } else {
      d=c(d,"black")
    }
    
    for (u in 2:h){
      c=toString(t[u,1])
      e=toString(t[(u-1),1])
      if (c != e) {
        if ((t[u,7]>1) | (t[u,3]=="hom")) {
          d=c(d,"red")
        } else {
          d=c(d,"black")
        }
      }
    }
    p2 <- p2 + theme(axis.text.y=element_text(colour=d))
    
    #size depending on variant occurence
    p2 <- p2 + scale_size_manual(values = c("1"= 5, "2"= 4, "3"= 4, "4"= 4, "5"= 4, "6"= 3, "7"= 3, "8"= 3, "9"= 3, "10"= 3, "11"= 2, "12"= 2, "13"= 2, "14"= 2, "15"= 2, "16"= 2, "17"= 2, "18"= 2, "19"= 2, "20"= 2, "21"= 2, "22"= 2, "23"= 2, "24"= 2, "25"= 2, "26"= 2, "27"= 2, "28"= 2, "29"= 2, "30"= 2, ">x"= 2))
    
    #p2 <- p2 + scale_size_manual(breaks=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",">x"),labels=c("1","[2-5]","[2-5]","[2-5]","[2-5]","[6-10]","[6-10]","[6-10]","[6-10]","[6-10]","[11-15]","[11-15]","[11-15]","[11-15]","[11-15]","[15-20]","[15-20]","[15-20]","[15-20]","[15-20]",">x"),values = c("1"= 5, "2"= 4, "3"= 4, "4"= 4, "5"= 4, "6"= 3, "7"= 3, "8"= 3, "9"= 3, "10"= 3, "11"= 2, "12"= 2, "13"= 2, "14"= 2, "15"= 2, "16"= 2, "17"= 2, "18"= 2, "19"= 2, "20"= 2, ">x"= 2))

    # on fait commencer l'axe x a 0
    p2 <- p2 + scale_x_continuous(expand = c(0, 0))
    # on enleve les "ticks" sur l'axe y
    p2 <- p2 + theme(axis.ticks.y = element_blank())
    p2
    outName <- paste("results/", g, ".png", sep="")
    ggsave(filename=outName, width=40,height=23, units="cm")
  }
}

graphIt(listGene)
