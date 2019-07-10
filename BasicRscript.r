######################################################################
######################################################################
####											
####			   MÉTODOS QUANTITATIVOS PARA				
####				O ESTUDO DE CORPORA				
####			     DE TEXTOS MULTILÍNGUES	
####
####					QUANTITATIVE METHODS FOR THE STUDY
####					OF CORPORA OF MULTILINGUAL TEXTS
####											
######################################################################
######################################################################

# Authors: Arthur de Melo Sá, Rodrigo Araújo e Castro

###################################################################
#											
#             			PRELIMINARY STEPS			
#											
###################################################################

# loads and installs packages

# openxlsx does NOT depend on Java, like other packages

if(!require("openxlsx")) { 		 # if this package is not installed
	install.packages("openxlsx") # install it
	require("openxlsx")			 # and load it
} else {					 # if it is loaded
	require("openxlsx")			 # load it
}

# same as p_load() from pacman package

require(pacman)
p_load(installr)


# updates R with updateR() function from installr package

require(installr)
updateR()


# removes all objects from memory
rm(list=ls())

# lists objects in memory
ls()

getwd()

# set directory to chosen folder
setwd(dir=choose.dir())

# example
setwd(dir="C:/Users/Rodrigo/Desktop")


require(xlsx)
library(xlsx)


#save.image(file=choose.files())
## Salvar o script pessoal dentro do mesmo diretório da imagem



###################################################################
#											
#                     		VECTORS			
#											
###################################################################


# creates vectors:
tradutores <- c("Paula","Roberto","Vanessa","Oséias","José")
idades <- c(35,54,19,62,43)
sexo <- c("F","M","F","M","M")
orações <- c(4567,1235,7894,4789,1459)
publicada <- c(T,T,F,T,F)


# accesses vector elements:
tradutores[4]
idades[4]
sexo[1:3]
orações[c(2,5)]
length(publicada)
publicada[length(publicada)]
tradutores[length(tradutores)]


# checks structure of vectors:
class(tradutores)
class(idades)
class(sexo)
sexo <- as.factor(sexo)
is.factor(sexo)
is.character(sexo)

class(publicada)

mode(tradutores)
mode(publicada)
mode(sexo)

length(idades)

str(idades)
str(sexo)

sort(idades)
sort(tradutores)

sort(idades,decreasing=TRUE)

#
#ls()
#

decresc = sort(idades,decreasing=TRUE)

##
rm(decresc)
#ls()
##

summary(idades)
summary(orações)
summary(sexo)
summary(publicada)


# basic operations:
mean(idades)
mean(orações)
mean(tradutores)

median(idades)
median(orações)
median(tradutores)

sum(orações)
sum(idades)

max(idades)
min(orações)

range(idades)
range(orações)

diff(range(idades))	##Amplitude
diff(range(orações))	##Amplitude

max(idades)-min(idades)
max(orações)-min(orações)

prod(idades)
prod(orações)


# adding elements to vectors:
tradutores1 <- c("Paula","Roberto","Vanessa","Oséias","José",88)
tradutores1

tradutores2 <- c(tradutores,88)
tradutores2

tradutores3 <- c(tradutores, c(77,88,99))
tradutores3

tradutores4 <- c(tradutores, c("Juca","Adonirão","Zulmira"))
tradutores4

tradutores5 <- c(tradutores,tradutores2)
tradutores5

# removes elements from vectors:
tradutores6 <- tradutores2[-6]
tradutores6

tradutores7 <- tradutores3[-c(6:8)]
tradutores7

###################################################################
#											
#             			MATRICES			
#											
###################################################################


# creates vectors:
dados1 <- c(seq(1,11,by=2))
dados1

dados2 <- c(seq(300,200,by=-20))
dados2


# combines vectors:
matriz1 <- cbind(dados1,dados2)
matriz1

matriz2 <- rbind(dados1,dados2)
matriz2

matriz3 <- t(matriz1)
matriz3


# accesses elements from matrices:
matriz1[3,]
matriz1[,2]

matriz1[3,1]
matriz1[4,3]
matriz1[3,2]

matriz1[2:5,]
matriz1[2:5,2]


# adds a column to the matrix:
dados3 <- c(rep(1:3,each=2))
dados3

matriz3 <- cbind(matriz1,dados3)
matriz3


# adds a row to the matrix:
dados4 <- c(5:7)
dados4

matriz4 <- rbind(matriz3,dados4)
matriz4

matriz4a <- rbind(matriz3,c(5:7))
matriz4a

# works with matrices:
class(matriz4)
is.matrix(matriz4)

mode(matriz4)

length(matriz4)
length(matriz4[,2])
length(matriz4[2,])

dim(matriz4)

str(matriz4)

summary(matriz4)
summary(t(matriz4))

mean(matriz4)
mean(matriz4[,3])
mean(matriz4[5,])

colMeans(matriz4)
rowMeans(matriz4)

colMeans(matriz4[,1:2])
rowMeans(matriz4[1:3,])
rowMeans(matriz4[,1:2])

median(matriz4)

mediana <- sort(matriz4)
mediana
length(mediana)
mediana[11]

median(matriz4[,1])
median(matriz4[6,])

sort(matriz4[,1])

sum(matriz4)
sum(matriz4[,2])
sum(matriz4[4,])

colSums(matriz4)
rowSums(matriz4)

colSums(matriz4[,2:3])
rowSums(matriz4[4:6,])

max(matriz4)
max(matriz4[,2])
max(matriz4[,1])

min(matriz4)
min(matriz4[,2])

range(matriz4)				## range
diff(range(matriz4))		## amplitude

max(matriz4)-min(matriz4)

range(matriz4[,2])			## range
diff(range(matriz4[,2]))	## amplitude

prod(matriz4)
prod(matriz4[3,])

which.max(matriz4)
matriz4[8]

which.min(matriz4)
matriz4[1]

which.max(matriz4[,1])
which.min(matriz4[5,])



###################################################################
#											
#                 		DATAFRAMES			
#											
###################################################################


# creates a dataframe:
tradutores <- c("Paula","Roberto","Vanessa","Oséias","José")

idades <- c(35,54,20,62,43)

sexo <- c("F","M","F","M","M")
class(sexo)
sexo <- as.factor(sexo)
is.factor(sexo)

orações <- c(4567,1235,7894,4789,1459)

publicada <- c(T,T,F,T,F)

info <- data.frame(tradutores,idades,sexo,orações,publicada,stringsAsFactors = F)
info

# check structure of dataframe:
class(info)
is.data.frame(info)
is.matrix(info)

mode(info)

class(info$idades)
class(info$sexo)
class(info$tradutores)

mode(info$publicada)
mode(info$orações)

length(info)

str(info)

summary(info)
summary(info$idades)
summary(info$publicada)


# accesses dataframe elements:
info$idades
info$tradutores

info$idades[4]
info$tradutores[1:3]


# adds vector (column) to dataframe:
data <- c("1984","2008",NA,"1995",NA)
data

info2 <- data.frame(info,data,stringsAsFactors = F)
info2

is.na(info2)
#is.factor()

info3 <- !is.na(info2)
info3


# removes vector from dataframe:
info4 <- info2[-5]
info4


# labels rows and columns:
colnames(info4) <- c("TRADUTOR","IDADE","SEXO","ORAÇÕES","PUBLICAÇÃO")
info4

names(info4) <- c("TRADUTOR","IDADE","SEXO","ORAÇÕES","PUBLICAÇÃO")
info4

rownames(info4) <- c("Paula","Roberto","Vanessa","Oséias","José")
info4

info5 <- info4[-1]
info5

rownames(info5)
colnames(info5)


# basic operations:
mean(info5$IDADE)
mean(info5$ORAÇÕES)
mean(info5$PUBLICAÇÃO)

info5$PUBLICAÇÃO <- as.factor(info5$PUBLICAÇÃO)
is.factor(info5$PUBLICAÇÃO)

str(info5$PUBLICAÇÃO)
summary(info5$PUBLICAÇÃO)

median(info5$IDADE)
median(info5$ORAÇÕES)

sum(info5$ORAÇÕES)

max(info5$IDADE)
min(info5$IDADE)

range(info5$IDADE)

diff(range(info5$IDADE))	##Amplitude

max(info5$IDADE)-min(info5$IDADE)

prod(info5$ORAÇÕES)

info5$ORAÇÕES[which.max(info5$IDADE)]
info5$ORAÇÕES[4]

info5$SEXO[which.min(info5$ORAÇÕES)]

rownames(info5)[which.min(info5$PUBLICAÇÃO)]

pos <- which.min(info5$PUBLICAÇÃO)
pos

rownames(info5)[pos]

info5$PUBLICAÇÃO[info5$IDADE > 43]
info5$ORAÇÕES[info5$IDADE < 43]
info5$SEXO[info5$IDADE <= 43]

rownames(info5)[info5$IDADE >= 43]
rownames(info5)[info5$IDADE == 43]
rownames(info5)[info5$IDADE != 43]



###################################################################
#											
#           DATAFRAMES, TABLES AND GRAPHS	
#											
###################################################################


# imports file as dataframe:
require(xlsx)

rm(list=ls())

ls()

getwd()

processos <- read.xlsx("bancobliss2.xlsx", 
				sheetIndex=1,
				startRow=1,
				colIndex=1:9,
				header=TRUE,
				stringsAsFactors=FALSE)
processos

class(processos)

str(processos)

mode(processos)

# visualizing dataframe:
processos[1:3,]
processos[,2:4]

processos[3,2:4]
processos[3,4]

processos["negative"]

processos$negative

processos[1,"negative"]

processos$negative[1:5]

head(processos)
head(processos,6)

tail(processos)
tail(processos,6)

processos2 <- edit(processos)	# doesn't save the changes
fix(processos)	# saves changes

dimnames(processos)

rownames(processos)
colnames(processos)

rownames(processos)[processos$traducao == "ORIGINAL"]
rownames(processos)[processos$traducao == "ESPANHOL"]

rownames(processos)[processos$traducao != "PORTUGUES"]


# checks dataframe:
is.vector(processos$traducao)

str(processos)

processos$traducao <- as.factor(processos$traducao)
is.factor(processos$traducao)

class(processos$traducao)

processos$TEXTO <- as.factor(processos$TEXTO)

str(processos)

mode(processos)

summary(processos)


# working with dataframes:

# creates column "publicado"
publicado = c(rep("SIM",5),rep("NAO",4),rep("SIM",2))

# casting (transforming) vector into factor
publicado = as.factor(publicado)

# creates column "1a_edicao"
prim_edicao = c(rep("NAO",2),rep("SIM",7),rep("NAO",2))

# updating object with logic column
processos = cbind(processos[,1:2],publicado,prim_edicao,processos[,3:9])

# analyses

processos[(processos$traducao == "ESPANHOL") &
		(processos$publicado == "NAO"),]

processos[(processos$traducao == "PORTUGUES") &
		(processos$publicado == "SIM"),]

processos$traducao == "ESPANHOL"

processos$publicado == "SIM"

processos$publicado  != "SIM"

processos[(processos$publicado != "SIM"),"TEXTO"]


# separates a dataset: FUNCTION split()
publicacao <- split(processos,processos$publicado)
publicacao

publicacao$NAO

mode(publicacao)

class(publicacao)

publicacao [[1]]
publicacao [[1]][1:5,]

tipo_publicacao <- publicacao[[1]]

tipo_publicacao

class(tipo_publicacao)


# selects a dataset: FUNCTION subset()
subconjunto1  <- subset(processos, traducao == "PORTUGUES", 
				select=c(traducao,publicado))
subconjunto1

subconjunto2 <- subset(processos, 
				traducao == "ESPANHOL" & publicado == "SIM", 
				select=c(traducao,publicado))

subconjunto2


# exports the dataset:
	# creates the folder "Resultados" in your directory
	dir.create("Resultados")

write.xlsx(publicacao[[1]],"Resultados\\Publicacao1.xlsx",col.names=TRUE,row.names=TRUE)

write.table(publicacao[[1]],"Resultados\\Publicacao2.txt",sep=",",col.names=TRUE,row.names=TRUE)

write.csv(publicacao[[1]],"Resultados\\Publicacao1.csv",col.names=TRUE,row.names=TRUE)

write.csv2(publicacao[[1]],"Resultados\\Publicacao2.csv",col.names=TRUE,row.names=TRUE)


# creates a contingency (counting) table with the data from ONE column of the dataframe:
tabela_publicado <- table(processos$publicado)
tabela_publicado

class(tabela_publicado)

mode(tabela_publicado)

str(tabela_publicado)

summary(tabela_publicado)


tabela_traducao <- table(processos$traducao)
tabela_traducao

str(tabela_traducao)

summary(tabela_traducao)

# creates a contingency (counting) table with the data from TWO column of the dataframe:
tabela_publicado_traducao <- table(processos$publicado,processos$traducao)
tabela_publicado_traducao

str(tabela_publicado_traducao)

summary(tabela_publicado_traducao)

tabela_publicado_prim_edicao <- table(processos$publicado,processos$prim_edicao)
tabela_publicado_prim_edicao

tabela_publicado_prim_edicao2 <- table(publicado=processos$publicado,prim_edicao=processos$prim_edicao)
tabela_publicado_prim_edicao2 

str(tabela_publicado_prim_edicao)

summary(tabela_publicado_prim_edicao)

# creates a contingency (counting) table with the data from THREE column of the dataframe:
tabela_TRES <- table(	processos$traducao,
				processos$publicado,
				processos$prim_edicao,
					dnn = c("Translation",
						"Published",
						"First Edition")
			  )
tabela_TRES

str(tabela_TRES)

summary(tabela_TRES)


tabela_TRESa <- ftable(processos$traducao,processos$publicado,processos$prim_edicao,
				dnn = c("Translation","Published","First Edition"))
tabela_TRESa

str(tabela_TRESa)

summary(tabela_TRESa)


# creates a frequency table:
tabela_frequencias1 <- prop.table(table(processos$traducao))
tabela_frequencias1


# creates a relative frequency table:
tabela_frequencias2 <- prop.table(table(processos$traducao,
						  processos$publicado),
							margin=1)
tabela_frequencias2


tabela_frequencias2b <- prop.table(table(processos$traducao,
						  processos$publicado),
							margin=2)
tabela_frequencias2b

colSums(tabela_frequencias2b)


tabela_frequencias3 <- prop.table(ftable(processos$traducao,
						   processos$publicado,
						   processos$prim_edicao,
				dnn = c("Translation",
					"Published",
					"First Edition"))	,margin=2)
tabela_frequencias3

colSums(tabela_frequencias3)


tabela_frequencias3b <- prop.table(ftable(processos$traducao,
						   processos$publicado,
						   processos$prim_edicao,
				dnn = c("Translation",
					"Published",
					"First Edition"))	,margin=1)
tabela_frequencias3b

colSums(tabela_frequencias3)



# creates graphs to represent the factors:

## dispersion:
plot(tabela_publicado)

plot(tabela_publicado,type="p")
plot(tabela_publicado,type="l")
plot(tabela_publicado,type="b")
plot(tabela_publicado,type="c")
plot(tabela_publicado,type="o")
plot(tabela_publicado,type="h")
plot(tabela_publicado,type="s")

plot(tabela_publicado, type="p", pch="#")

plot(tabela_publicado,type="p",
	col=c("green","red","blue","yellow","purple","magenta"))

plot(tabela_publicado,type="h",
	col=c("green","red","blue","yellow","purple","magenta"))

colors()

plot(tabela_publicado,type="p",
	col=c("green","red","blue","yellow","purple","magenta"),
	main="GRÁFICO DE DISPERSÃO",
	sub="PUBLICACAO",
	xlab="Publicacao",
	ylab="Ocorrências"
	)

# export the graphs!!!

## Pie:
pie(tabela_publicado)

dev.new()
pie(tabela_traducao)
dev.off()

dev.new()

tabela_prim_edicao = table(processos$prim_edicao)
pie(tabela_prim_edicao, main="gráfico de setores (pizza)")

dev.new()
par(mfrow=c(2,2))
pie(tabela_prim_edicao)
pie(tabela_prim_edicao)
pie(tabela_prim_edicao)

par(mfrow=c(1,1))
dev.off()

# Export the graphs!!!


## Bars:
	# This function does NOT accept factors as arguments
	# we can only give a "table" object
dev.new()
barplot(tabela_publicado)

dev.new()
barplot(tabela_publicado, col=c("green","red"))

dev.new()
barplot(tabela_traducao,
		main="Gráfico de barras",
		col=c("green","red","yellow","blue"),
		names=c("GM","VMSC","VSD"),
		xlab="Tipo de verbo",
		ylab="Número de ocorrências",
		cex.names=0.75,					# column size
		cex.lab=2,						# axes sizes
		)

dev.new()
barplot(tabela_publicado,
		legend=rownames(tabela_publicado),
		col=c("green","red","yellow","blue"),
		beside=TRUE,
		offset=TRUE
		)

dev.new()
barplot(tabela_prim_edicao,
		legend=rownames(tabela_publicado),
		col=c("green","red","yellow","blue"),
		beside=TRUE,
		offset=TRUE,
		horiz=TRUE
		)

# export graphs!!!


#########################################
# imports a new file as dataframe:
# (not included with the script, belongs to the LETRA LAB FROM UFMG)
bliss <- read.csv("bancobliss.csv", sep=";")
bliss

class(bliss)

mode(bliss)

str(bliss)


install.packages("psych")
require(psych)


# selects the columns "Words" and "Sentences",
#	grouped by linguistic system: FUNCTION split()
words_system <- split(bliss$words,bliss$traducao)
words_system

class(words_system)

mode(words_system)

sentences_system <- split(bliss$sentences,bliss$traducao)
sentences_system

class(sentences_system)

mode(sentences_system)


# describes the variables "Words" and "Sentences": FUNCTION describeBy()
words <- describeBy(bliss$words,bliss$traducao)   ; words
sentences <- describeBy(bliss$sentence,bliss$traducao)   ; sentences


# summarizes the variables "Words" and "Sentences"
#	according to linguistic system: FUNCTION aggregate()
word.sum <- aggregate(bliss$words,
		by=list(bliss$traducao),
            FUN=summary)
word.sum 


sentence.sum <- aggregate(bliss$sentence,
			by=list(bliss$traducao),
			FUN=summary)
sentence.sum


# creates a dataframe with the mean and the standard deviation of the variables
#	"Words" and "Sentences" for the translation into PORTUGUESE (portugues):
portugues <- words$PORTUGUES[3:4] 
portugues


portugues <- rbind(portugues,sentences$PORTUGUES[3:4])
portugues 


portugues <- data.frame(c("words","sentences"),portugues)
portugues
names(portugues) <- c("Variável","Media","Desvio-Padrão")
portugues


class(portugues)

mode(portugues)


# creates a dataframe with the mean and the standard deviation of the variables
#	"Words" and "Sentences" for the translation into SPANISH (espanhol):
espanhol <- words$ESPANHOL[3:4]

espanhol <- rbind(espanhol,sentences$ESPANHOL[3:4])

espanhol <- data.frame(c("words","sentences"),espanhol)

names(espanhol) <- names(portugues)

espanhol


# compares the translations through the error bars: FUNCTION error.bars()
bliss_portugues <- subset(bliss, bliss$traducao=="PORTUGUES")
bliss_portugues

bliss_espanhol <- subset(bliss, bliss$traducao=="ESPANHOL")
bliss_espanhol

words_bliss <- data.frame(bliss_portugues$words,bliss.espanhol$words)
words_bliss

names(words_bliss) <- c("Portugues","Espanhol")
words_bliss

class(words_bliss)

# creating the graph:
error.bars(words_bliss,
          xlab="Tradução", ylab="Número de palavras", bar=TRUE)

# customizing the graph:
error.bars(words_bliss, xlab="Tradução", ylab="Número de palavras", 
          bar=TRUE, main="Intervalos de 95% de confiança", arrow.len=0.2)


# create general graphs:

## dispersion:
plot(words_bliss$Portugues)

plot(words_bliss$Portugues,type="p")
plot(words_bliss$Portugues,type="l")
plot(words_bliss$Portugues,type="b")
plot(words_bliss$Portugues,type="c")
plot(words_bliss$Portugues,type="o")
plot(words_bliss$Portugues,type="h")
plot(words_bliss$Portugues,type="s")

plot(words_bliss$Portugues, type="p", pch="#")

plot(words_bliss$Portugues,type="p",
	col=c("green","red","blue","purple","magenta"))
plot(words_bliss$Portugues,type="h",
	col=c("green","red","blue","purple","magenta"))

plot(words_bliss$Portugues,type="p",
	col=c("green","red","blue","purple","magenta"),
	main="GRÁFICO DE DISPERSÃO",
	sub="NÚMERO DE PALAVRA",
	xlab="Tradução",
	ylab="Número de palavras"
	)

# export the graphs!!!


## Pie:
positive_negative_portugues <- cbind(bliss$positive[7:11],bliss$negative[7:11])
colnames(positive_negative_portugues) <- c("positive","negative")
positive_negative_portugues

class(positive_negative_portugues)

mode(positive_negative_portugues)


media_pos_neg <- colMeans(positive_negative_portugues)
media_pos_neg

proportions <- round(prop.table(media_pos_neg),2)
proportions

dev.new()
pie(media_pos_neg, main="gráfico de setores (pizza)",
		labels=c("Positive - 87%","Negative - 13%"),
		col=c("yellow","white"),
		border=NA,
		init.angle=20)

# export the graphs!!!


## Bars:
dev.new()
barplot(positive_negative_portugues)

dev.new()
barplot(positive_negative_portugues, col=c("green","red","yellow","blue","purple"))

dev.new()
barplot(positive_negative_portugues,
		beside=TRUE,
		main="Gráfico de barras",
		col=c("green","red","yellow","blue","purple"),
		legend=c("ACC","EV","EVS","JC","MS"),
		names=c("Positive","Negative"),
		xlab="Tipo de verbo",
		ylab="Número de ocorrências"
		)

# export the graphs!!!



###################################################################
#											
#           	  	FUNCTIONS		
#											
###################################################################


# functions that we have seen so far:
mean()
median()
min()
max()
diff(range())


# creates a function to perform these operations:

my_function <- function(objeto)
{
	my_mean	 <-	mean(objeto)
	my_median <- median(objeto)
	my_min	 <-	min(objeto)
	my_max	 <-	max(objeto)
	my_amplitude <-	diff(range(objeto))
	out		 <-	c(my_mean, my_median, my_min, my_max, my_amplitude)
	names(out)	 <-	c("Média", "Mediana", "Mínimo", "Máximo", "Amplitude")
	out
}


# applies the functions on any object:
# example = sampling of numbers
x <- sample(89:4631)

my_function(x)


# Opening a new script:
#	(it is possible to open several scripts in R screen)
#
#		1) Copy and paste the function above in a new script;
#		2) Save the script in your directory with the name "my_script.R"
#		3) Close the script that you just saved.


# opens a new script in the background:
source("my_script.R")

ls()


# checks if the created function is working:
my_function

my_function(x)

#######################################################
