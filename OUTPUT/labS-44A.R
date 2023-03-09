##########################

LABORATORIO 

###########################
ADICIONALES
library(dplyr)
library(ggplot2)
library(stats)


#CARGAR LIBRERIA NECESARIA
library(stats)


#Attaching package: 'gplots'
#lowess

#ANEXAMOS DATOS 
h1 <- c(10,20,10,20,10,20,10,20)
h2 <- c(20,10,20,10,20,10,20,10)

l1 <- c(1,3,1,3,1,3,1,3)
l2 <- c(3,1,3,1,3,1,3,1)

mat <- rbind(h1,h2,l1,l2)

#CREAR GRÁFICO

par(mfrow =c(1,1), mar=c(4,4,1,1))
plot(1:8,rep(0,8), ylim=c(0,35), pch="", xlab="Time", ylab="Gene Expression")
for (i in 1:nrow(mat)) {lines(1:8,mat[i,], lwd=3, col=i)}
legend(1,35,rownames(mat), 1:4, cex=0.7)

#CALCULAR DISTANCIA 
dist(mat)

#SACAR PLOT
## I will use the default for linkage method: complete
plot(hclust(dist(mat)))

#CALCULAR HEATMAP DE RENGLONES
library(ggplot2)

heatmap (mat, Colv=NA, scale = "row")
heatmap (mat, Colv=NA, scale = "none")

mat.scaled<- t(scale(t(mat), center=TRUE, scale = TRUE))
mat.scaled

#DISTANCIA ENTRE GENES

dist(mat.scaled)

#NUEVO GRAFICO
plot(hclust(dist(mat.scaled)))

heatmap(mat.scaled, Colv = NA, scale = "none")

#CORRELACION ENTRE GENES
cor(t(mat))

#DISTANCIA EN LAS CORRELACIONES
1- cor(t(mat))

#GRAFICAS
hc <- hclust(as.dist(1-cor(t(mat))))
plot(hc)

heatmap(mat, Colv = NA, Rowv=as.dendrogram(hc), scale = "none")

#SEGUNDA oPCION DE HEATMAPS
heatmap(mat, trace = "none", Colv= NA, dendrogram = "row", scale = "none")

#TERCERA OPCION HEATMAPS
heatmap(mat, trace = "none", Colv= NA, dendrogram = "row", scale = "row")
        
#CUARTA OPCION HEATMAPS USANDO DISTANCIA ECLUDIANA
heatmap (t(scale(t(mat), center=TRUE, scale=TRUE)), trace = "none", Colv= NA, dendrogram = "row")

#USO 1-COR (x) DE CORRELACION EN HEAT
heatmap (mat, trace = "none", 
  Colv= NA, dendrogram = "row",
  scale = "none",
  hclust=function(x) hclust(x, method="complete"), distfun=function(x) as.dist(1-cor(t(x))))

#OTRA ESCALA REPRESENTADA POR LOS COLORES
heatmap(mat, trace = "none", 
          Colv= NA, dendrogram = "row",
          scale = "row",
          hclust=function(x) hclust(x, method='complete'), distfun=function(x) as.dist(1-cor(t(x))))                                                                              

#ULTIMA COLORACION DE CLUSTERS
heatmap(t(scale(t(mat), center=TRUE, scale=TRUE)), trace = "none", 
Colv= NA, dendrogram = "row",
hclust=function(x) hclust(x, method='complete'), distfun=function(x) as.dist(1-cor(t(x))))
