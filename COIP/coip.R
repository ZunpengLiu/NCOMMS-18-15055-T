library(readr)
DGCR8_CHAPS_NP40_1 <- read_csv("/Volumes/Zunpeng/08_Result/04_DGCR8/质谱/DGCR8-CHAPS-NP40_1.csv")
View(DGCR8_CHAPS_NP40_1)

res<-DGCR8_CHAPS_NP40_1
res<-res[-1,]
res$lgscore<-log2(res$Score)

pdf("Liver malonyllysine point plot.pdf", width = 8, height = 6)
ggplot(res,aes(y=lgscore,x=res$Coverage,
                  #colour="Subcellular localization"，
                  color=res$color,size=res$Score))+
  geom_point()+
  scale_color_manual(values=c("dodgerblue3","grey","brown2"))+
  theme_few() +  
  labs(x="Coverage (%)",y="log2(score)")+
  geom_label_repel(aes(res$Coverage, res$lgscore,label=res$label,colour=factor(res$color)),label.size = 0.05,
                   fill="white")
dev.off()











pdf("Liver malonyllysine point plot.pdf", width = 8, height = 6)
ggplot(res,aes(y=lgscore,x=sort,
               #colour="Subcellular localization"，
               color=res$color,size=res$Score))+
  geom_point()+
  scale_color_manual(values=c("dodgerblue3","grey","brown2"))+
  theme_few() +  
  labs(x="sort",y="log2(score)")+
  geom_label_repel(aes(res$sort, res$lgscore,label=res$label,colour=factor(res$color)),label.size = 0.05,
                   fill="white")
dev.off()














