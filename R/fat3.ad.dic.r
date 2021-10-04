#' Fatorial triplo com um tratamento adicional em DIC
#'
#' \code{fat3.ad.dic} Analisa experimentos em fatorial triplo
#' com um tratamento adicional em Delineamento Inteiramente
#' Casualizado balanceado, considerando o modelo fixo.
#' @param fator1 Vetor numerico ou complexo contendo os niveis
#' do fator 1.
#' @param fator2 Vetor numerico ou complexo contendo os niveis
#' do fator 2.
#' @param fator3 Vetor numerico ou complexo contendo os niveis
#' do fator 3.
#' @param repet Vetor numerico ou complexo contendo as
#' repeticoes.
#' @param resp Vetor numerico ou complexo contendo a variavel
#' resposta.
#' @param respAd Vetor numerico ou complexo contendo a variavel
#' resposta do tratamento adicional.
#' @param quali Logico, se TRUE (default) na primeira posicao,
#' os niveis do fator 1 sao entendidos como qualitativos, se
#' FALSE, quantitativos; da mesma forma, a segunda posicao e
#' referente aos niveis do fator 2, e a terceira aos niveis do
#' fator 3.
#' @param mcomp Permite escolher o teste de comparacao multipla;
#' o \emph{default} e o teste de Tukey, contudo tem-se como
#' outras opcoes: o teste LSD ('lsd'), o teste LSDB ('lsdb'),
#' o teste de Duncan ('duncan'), o teste de SNK ('snk'), o
#' teste de Scott-Knott ('sk'), o teste de comparacoes
#' multiplas bootstrap ('ccboot') e o teste de Calinski e
#' Corsten baseado na distribuicao F ('ccf').
#' @param fac.names Permite nomear os fatores 1, 2 e 3.
#' @param sigT Significancia a ser adotada pelo teste de
#' comparacao multipla de medias; o default e 5\%.
#' @param sigF Significancia a ser adotada pelo teste F da
#' ANAVA; o default e 5\%.
#' @param unfold Orienta os desdobramentos apos a analise de
#' variancia. Se NULL (\emph{default}), sao feitas as analises
#' recomendadas; se '0', e feita apenas a analise de variancia;
#' se '1', os efeitos simples sao estudados; se '2.1', '2.2' ou
#' '2.3, as interaoes duplas 1, 2 ou 3 sao estudadas; se '3',
#' a interacao tripla e estudada.
#' @details Os argumentos sigT e mcomp so serao utilizados
#' quando os tratamentos forem qualitativos.
#' @return Sao retornados os valores da analise de variancia
#' do DIC em questao com um tratamento adicional,o teste de
#' normalidade de Shapiro-Wilk para os residuos do modelo, o
#' ajuste de modelos de regressao (caso de tratamentos
#' quantitativos) ou os testes de comparacao de medias (caso de
#' tratamentos qualitativos): teste de Tukey, teste de Duncan,
#' teste t de Student (LSD), teste t de Bonferroni, teste de
#' Student-Newman-Keuls (SNK), teste de Scott-Knott e teste de
#' comparacoes multiplas bootstrap; com o desdobramento da
#' interacao, caso esta seja significativa.
#' @references HEALY, M. J. R. The analysis of a factorial
#' experiment with additional treatments. Journal of
#' Agricultural Science, Cambridge, v. 47, p. 205-206.
#' 1956.
#' @author Eric B Ferreira,
#'\email{eric.ferreira@@unifal-mg.edu.br}
#' @author Denismar Alves Nogueira
#' @author Portya Piscitelli Cavalcanti
#' @note O \code{\link{graficos}} pode ser usado para
#' construir os graficos da regressao e o
#' \code{\link{plotres}} para analise do residuo da anava.
#' @seealso \code{\link{fat2.dic}}, \code{\link{fat2.dbc}},
#' \code{\link{fat3.dic}}, \code{\link{fat3.dbc}},
#' \code{\link{fat2.ad.dic}}, \code{\link{fat2.ad.dbc}},
#' and \code{\link{fat3.ad.dbc}}.
#' @examples
#' data(ex6)
#' attach(ex6)
#' data(respAd)
#' fat3.ad.dic(fatorA, fatorB, fatorC, rep, resp, respAd,
#' quali = c(TRUE, TRUE, TRUE), mcomp = "duncan", fac.names =
#' c("Fator A", "Fator B", "Fator C"), sigT=0.05, sigF = 0.05,
#' unfold=NULL)
#' @export

fat3.ad.dic<-function(fator1,
fator2,
fator3,
repet,
resp,
respAd,
quali=c(TRUE,TRUE,TRUE),
mcomp='tukey',
fac.names=c('F1','F2','F3'),
sigT=0.05,
sigF=0.05,
unfold=NULL) {

cat('------------------------------------------------------------------------\nLegenda:\n')
cat('FATOR 1: ',fac.names[1],'\n')
cat('FATOR 2: ',fac.names[2],'\n')
cat('FATOR 3: ',fac.names[3],'\n------------------------------------------------------------------------\n\n')

fatores<-data.frame(fator1,fator2,fator3)
Fator1<-factor(fator1)
Fator2<-factor(fator2)
Fator3<-factor(fator3)
nv1<-length(summary(Fator1)) #Diz quantos niveis tem o fator 1.
nv2<-length(summary(Fator2)) #Diz quantos niveis tem o fator 2.
nv3<-length(summary(Fator3)) #Diz quantos niveis tem o fator 3.
J<-(length(resp))/(nv1*nv2*nv3)
n.trat2<-nv1*nv2
n.trat3<-nv1*nv2*nv3
lf1<-levels(Fator1)
lf2<-levels(Fator2)
lf3<-levels(Fator3)

#Anava do fatorial 3
anavaF3<-summary(aov(resp~Fator1*Fator2*Fator3))
SQa<-anavaF3[[1]][1,2]
SQb<-anavaF3[[1]][2,2]
SQc<-anavaF3[[1]][3,2]
SQab<-anavaF3[[1]][4,2]
SQac<-anavaF3[[1]][5,2]
SQbc<-anavaF3[[1]][6,2]
SQabc<-anavaF3[[1]][7,2]
gla=nv1-1
glb=nv2-1
glc=nv3-1
glab=(nv1-1)*(nv2-1)
glac=(nv1-1)*(nv3-1)
glbc=(nv2-1)*(nv3-1)
glabc=(nv1-1)*(nv2-1)*(nv3-1)
QMa=SQa/gla
QMb=SQb/glb
QMc=SQc/glc
QMab=SQab/glab
QMac=SQac/glac
QMbc=SQbc/glbc
QMabc=SQabc/glabc

#Anava de todos os tratamentos do experimento (fatorial 3 + adicional)
col1<-numeric(0)
for(i in 1:n.trat3) {
col1<-c(col1, rep(i,J))
 }
col1<-c(col1,rep('ad',J))
col2<-c(repet,rep(1:J))
col3<-c(resp,respAd)
tabF3ad<-data.frame("TRAT"=col1, "REP"=col2, "RESP2"=col3)
TRAT<-factor(tabF3ad[,1])
anava<-aov(tabF3ad[,3]~TRAT)
anavaTr<-summary(anava)
SQad<-anavaTr[[1]][1,2] - (SQa+SQb+SQc+SQab+SQac+SQbc+SQabc)
SQE<-anavaTr[[1]][2,2]
SQT<-anavaTr[[1]][1,2]+anavaTr[[1]][2,2]
glad=1
glT=(nv1*nv2*nv3+1)*J-1
glE=glT-(gla+glb+glc+glab+glac+glbc+glabc+1)
QMad=SQad/glad
QME=SQE/glE
QMT=SQT/glT
Fca=QMa/QME
Fcb=QMb/QME
Fcc=QMc/QME
Fcab=QMab/QME
Fcac=QMac/QME
Fcbc=QMbc/QME
Fcabc=QMabc/QME
Fcad=QMad/QME
pv.fs=c(1-pf(Fca,gla,glE), 1-pf(Fcb,glb,glE), 1-pf(Fcc,glc,glE))

#Montando a tabela da ANAVA
an<-data.frame("GL"=c(gla, glb, glc, glab, glac, glbc, glabc, glad, glE, glT ),
"SQ"=c(round(c(SQa,SQb,SQc,SQab,SQac,SQbc,SQabc,SQad,SQE,SQT),5)),
"QM"=c(round(c(QMa,QMb,QMc,QMab,QMac,QMbc,QMabc,QMad,QME),5),''),
"Fc"=c(round(c(Fca,Fcb,Fcc,Fcab,Fcac,Fcbc,Fcabc,Fcad),4),'',''),
"Pr>Fc"=c(round(c(pv.fs, 1-pf(Fcab,glab,glE), 1-pf(Fcac,glac,glE), 1-pf(Fcbc,glbc,glE), 1-pf(Fcabc,glabc,glE), 1-pf(Fcad,glad,glE)),4), ' ', ' '))
colnames(an)[5]="Pr>Fc"
rownames(an)=c(fac.names[1],fac.names[2],fac.names[3],paste(fac.names[1],'*',fac.names[2],sep=''),paste(fac.names[1],'*',fac.names[3],sep=''),
paste(fac.names[2],'*',fac.names[3],sep=''),paste(fac.names[1],'*',fac.names[2],'*',fac.names[3],sep=''),"Ad vs Fatorial","Residuo","Total")
cat('------------------------------------------------------------------------
Quadro da analise de variancia\n------------------------------------------------------------------------\n')
print(an)
cat('------------------------------------------------------------------------\n')
pvalor<-c(1-pf(Fca,gla,glE), 1-pf(Fcb,glb,glE), 1-pf(Fcc,glc,glE), 1-pf(Fcab,glab,glE), 1-pf(Fcac,glac,glE), 1-pf(Fcbc,glbc,glE), 1-pf(Fcabc,glabc,glE))
#CV
cv<-round(sqrt(QME)/mean(col3)*100, 2)
cat('CV =',cv,'%\n')

#Teste de normalidade
pvalor.shapiro<-shapiro.test(anava$residuals)$p.value
cat('\n------------------------------------------------------------------------
Teste de normalidade dos residuos (Shapiro-Wilk)\n')
cat('valor-p: ',pvalor.shapiro, '\n')
if(pvalor.shapiro<0.05){cat('ATENCAO: a 5% de significancia, os residuos nao podem ser considerados normais!
------------------------------------------------------------------------\n')}
else{cat('De acordo com o teste de Shapiro-Wilk a 5% de significancia, os residuos podem ser considerados normais.
------------------------------------------------------------------------\n\n')}

#Contraste Ad vs Fatorial
cat('Contraste do tratamento adicional com o fatorial
------------------------------------------------------------------------\n')
x<-mean(respAd)
y<-mean(resp)

if(1-pf(Fcad,glad,glE)>sigF) { C1<-data.frame("Medias"=c(x,y))
rownames(C1)=c("Adicional","Fatorial")
colnames(C1)<-c("Medias")
cat('De acordo com o teste F, as medias dos dois grupos sao estatisticamente iguais.\n')
print(C1) }else{
C2<-data.frame("Media"=c(x,y),
" "=c(letters[1],letters[2]))
rownames(C2)=c("Adicional","Fatorial")
colnames(C2)<-c("Medias"," ")
print(C2)
}
cat('------------------------------------------------------------------------\n')

# Creating unfold #########################################
if(is.null(unfold)){
if(1-pf(Fcab,glab,glE)>sigF &&
 1-pf(Fcac,glac,glE)>sigF &&
 1-pf(Fcbc,glbc,glE)>sigF &&
 1-pf(Fcabc,glabc,glE)>sigF) {unfold<-c(unfold,1)}
if(1-pf(Fcabc,glabc,glE)>sigF &&
 1-pf(Fcab,glab,glE)<=sigF){unfold<-c(unfold,2.1)}
if(1-pf(Fcabc,glabc,glE)>sigF &&
 1-pf(Fcac,glac,glE)<=sigF){unfold<-c(unfold,2.2)}
if(1-pf(Fcabc,glabc,glE)>sigF &&
 1-pf(Fcbc,glbc,glE)<=sigF){unfold<-c(unfold,2.3)}
if(1-pf(Fcabc,glabc,glE)<=sigF){unfold<-c(unfold,3)}
}

#Para nenhuma interacao significativa, fazer...
if(any(unfold==1)) {
cat('\nInteracao nao significativa: analisando os efeitos simples
------------------------------------------------------------------------\n')
fatores<-data.frame('fator 1'=fator1,'fator 2' = fator2,'fator 3' = fator3)

for(i in 1:3){
#Para os fatores QUALITATIVOS, teste de Tukey
if(quali[i]==TRUE && pvalor[i]<=sigF) {
cat(fac.names[i])
if(mcomp=='tukey'){tukey(resp,fatores[,i],an[9,1],an[9,2],sigT)}
if(mcomp=='duncan'){duncan(resp,fatores[,i],an[9,1],an[9,2],sigT)}
if(mcomp=='lsd'){lsd(resp,fatores[,i],an[9,1],an[9,2],sigT)}
if(mcomp=='lsdb'){lsdb(resp,fatores[,i],an[9,1],an[9,2],sigT)}
if(mcomp=='sk'){scottknott(resp,fatores[,i],an[9,1],an[9,2],sigT)}
if(mcomp=='snk'){snk(resp,fatores[,i],an[9,1],an[9,2],sigT)}
if(mcomp=="ccboot"){ccboot(resp,fatores[,i],an[9,1],an[9,2],sigT)}
if(mcomp=="ccF"){ccF(resp,fatores[,i],an[9,1],an[9,2],sigT)}
}
if(quali[i]==TRUE && pvalor[i]>sigF) {
cat(fac.names[i])
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp,fatores[,i],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------')
}

#Para os fatores QUANTITATIVOS, regressao
if(quali[i]==FALSE && pvalor[i]<=sigF){
cat(fac.names[i])
reg.poly(resp, fatores[,i], an[9,1], an[9,2], an[i,1], an[i,2])
}

if(quali[i]==FALSE && pvalor[i]>sigF) {
cat(fac.names[i])
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp,fatores[,i],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------')
}
cat('\n')
}
}

#Se a(s) interacao(oes) dupla(s) for(em) significativa(s), desdobramento:

#Interacao Fator1*Fator2
if(any(unfold==2.1)){
cat("\n\n\nInteracao",paste(fac.names[1],'*',fac.names[2],sep='')," significativa: desdobrando a interacao
------------------------------------------------------------------------\n")

#Desdobramento de FATOR 1 dentro do niveis de FATOR 2
cat("\nDesdobrando ", fac.names[1], ' dentro de cada nivel de ', fac.names[2], '
------------------------------------------------------------------------\n')

des1<-aov(resp~Fator2/Fator1)

l1<-vector('list',nv2)
names(l1)<-names(summary(Fator2))
v<-numeric(0)
for(j in 1:nv2) {
for(i in 0:(nv1-2)) v<-cbind(v,i*nv2+j)
l1[[j]]<-v
v<-numeric(0)
}
des1.tab<-summary(des1,split=list('Fator2:Fator1'=l1))[[1]]

#Montando a tabela de ANAVA do des1
glf1=c(as.numeric(des1.tab[3:(nv2+2),1]))

SQf1=c(as.numeric(des1.tab[3:(nv2+2),2]))

QMf1=SQf1/glf1

Fcf1=QMf1/QME

rn<-numeric(0)
for(i in 1:nv2){ rn<-c(rn, paste(paste(fac.names[1],':',fac.names[2],sep=''),lf2[i]))}

anavad1<-data.frame("GL"=c(glf1, glE),
"SQ"=c(round(c(SQf1,SQE),5)),
"QM"=c(round(c(QMf1,QME),5)),
"Fc"=c(round(Fcf1,4),''),
"Pr>Fc"=c(round(1-pf(Fcf1,glf1,glE),4),' '))
colnames(anavad1)[5]="Pr>Fc"
rownames(anavad1)=c(rn,"Residuo")
cat('------------------------------------------------------------------------
Quadro da analise de variancia\n------------------------------------------------------------------------\n')
print(anavad1)
cat('------------------------------------------------------------------------\n\n')

ii<-0
for(i in 1:nv2) {
ii<-ii+1
if(1-pf(Fcf1,glf1,glE)[ii]<=sigF){
if(quali[1]==TRUE){
cat('\n\n',fac.names[1],' dentro do nivel ',lf2[i],' de ',fac.names[2],'
------------------------------------------------------------------------')
if(mcomp=='tukey'){
tukey(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=='duncan'){
duncan(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='lsd'){
lsd(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=='lsdb'){
lsdb(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='sk'){
scottknott(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='snk'){
snk(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=="ccboot"){
ccboot(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=="ccF"){
ccF(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
 }
}
else{#regressao
cat('\n\n',fac.names[1],' dentro do nivel ',lf2[i],' de ',fac.names[2],'
------------------------------------------------------------------------')
reg.poly(resp[Fator2==lf2[i]], fator1[Fator2==lf2[i]], an[9,1], an[9,2], des1.tab[i+2,1], des1.tab[i+2,2])
}
}
else{cat('\n\n',fac.names[1],' dentro do nivel ',lf2[i],' de ',fac.names[2],'\n')
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------\n')
}
 }
cat('\n\n')

#Desdobramento de FATOR 2 dentro do niveis de FATOR 1
cat("\nDesdobrando ", fac.names[2], ' dentro de cada nivel de ', fac.names[1], '
------------------------------------------------------------------------\n')

des2<-aov(resp~Fator1/Fator2)

l2<-vector('list',nv1)
names(l2)<-names(summary(Fator1))
v<-numeric(0)
for(j in 1:nv1) {
for(i in 0:(nv2-2)) v<-cbind(v,i*nv1+j)
l2[[j]]<-v
v<-numeric(0)
}
des2.tab<-summary(des2,split=list('Fator1:Fator2'=l2))[[1]]

#Montando a tabela de ANAVA do des2
glf2=c(as.numeric(des2.tab[3:(nv1+2),1]))

SQf2=c(as.numeric(des2.tab[3:(nv1+2),2]))

QMf2=SQf2/glf2

Fcf2=QMf2/QME

rn<-numeric(0)
for(k in 1:nv1){ rn<-c(rn, paste(paste(fac.names[2],':',fac.names[1],sep=''),lf1[k]))}

anavad2<-data.frame("GL"=c(glf2, glE),
"SQ"=c(round(c(SQf2,SQE),5)),
"QM"=c(round(c(QMf2,QME),5)),
"Fc"=c(round(Fcf2,4),''),
"Pr>Fc"=c(round(1-pf(Fcf2,glf2,glE),4),' '))
colnames(anavad2)[5]="Pr>Fc"
rownames(anavad2)=c(rn,"Residuo")
cat('------------------------------------------------------------------------
Quadro da analise de variancia\n------------------------------------------------------------------------\n')
print(anavad2)
cat('------------------------------------------------------------------------\n\n')

ii<-0
for(i in 1:nv1) {
ii<-ii+1
if(1-pf(Fcf2,glf2,glE)[ii]<=sigF){
if(quali[2]==TRUE){
cat('\n\n',fac.names[2],' dentro do nivel ',lf1[i],' de ',fac.names[1],'
------------------------------------------------------------------------')
if(mcomp=='tukey'){
tukey(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=='duncan'){
duncan(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='lsd'){
lsd(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=='lsdb'){
lsdb(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='sk'){
scottknott(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='snk'){
snk(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=="ccboot"){
ccboot(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=="ccF"){
ccF(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
 }
}
else{#regressao
cat('\n\n',fac.names[2],' dentro do nivel ',lf1[i],' de ',fac.names[1],'
------------------------------------------------------------------------')
reg.poly(resp[Fator1==lf1[i]], fator2[Fator1==lf1[i]], an[9,1], an[9,2], des2.tab[i+2,1], des2.tab[i+2,2])
}
}
else{cat('\n\n',fac.names[2],' dentro do nivel ',lf1[i],' de ',fac.names[1],'\n')
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------\n')
}
}

#Checar o Fator3
if(pvalor[5]>sigF && pvalor[6]>sigF) {
cat('\nAnalisando os efeitos simples do fator ',fac.names[3],'
------------------------------------------------------------------------\n')

i<-3
{
#Para os fatores QUALITATIVOS, teste de Tukey
if(quali[i]==TRUE && pvalor[i]<=sigF) {
 cat(fac.names[i])
if(mcomp=='tukey'){tukey(resp,fatores[,i],an[8,1],an[8,2],sigT)}
if(mcomp=='duncan'){duncan(resp,fatores[,i],an[8,1],an[8,2],sigT)}
if(mcomp=='lsd'){lsd(resp,fatores[,i],an[8,1],an[8,2],sigT)}
if(mcomp=='lsdb'){lsdb(resp,fatores[,i],an[8,1],an[8,2],sigT)}
if(mcomp=='sk'){scottknott(resp,fatores[,i],an[8,1],an[8,2],sigT)}
if(mcomp=='snk'){snk(resp,fatores[,i],an[8,1],an[8,2],sigT)}
if(mcomp=="ccboot"){ccboot(resp,fatores[,i],an[8,1],an[8,2],sigT)}
if(mcomp=="ccF"){ccF(resp,fatores[,i],an[8,1],an[8,2],sigT)}
}

if(quali[i]==TRUE && pvalor[i]>sigF) {
cat(fac.names[i])
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp,fatores[,i],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------')
}

#Para os fatores QUANTITATIVOS, regressao
if(quali[i]==FALSE && pvalor[i]<=sigF){
cat(fac.names[i])
reg.poly(resp, fatores[,i], an[8,1],an[8,2], an[i,1], an[i,2])
}

if(quali[i]==FALSE && pvalor[i]>sigF) {
cat(fac.names[i])
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp,fatores[,i],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------')
}
cat('\n')
}
}
}

#Interacao Fator1*Fator3
if(any(unfold==2.2)){
cat("\n\n\nInteracao",paste(fac.names[1],'*',fac.names[3],sep='')," significativa: desdobrando a interacao
------------------------------------------------------------------------\n")

#Desdobramento de FATOR 1 dentro do niveis de FATOR 3
cat("\nDesdobrando ", fac.names[1], ' dentro de cada nivel de ', fac.names[3], '
------------------------------------------------------------------------\n')

des3<-aov(resp~Fator3/Fator1)

l1<-vector('list',nv3)
names(l1)<-names(summary(Fator3))
v<-numeric(0)
for(j in 1:nv3) {
for(i in 0:(nv1-2)) v<-cbind(v,i*nv3+j)
l1[[j]]<-v
v<-numeric(0)
}
des3.tab<-summary(des3,split=list('Fator3:Fator1'=l1))[[1]]

#Montando a tabela de ANAVA do des3
glf3=c(as.numeric(des3.tab[3:(nv3+2),1]))

SQf3=c(as.numeric(des3.tab[3:(nv3+2),2]))

QMf3=SQf3/glf3

Fcf3=QMf3/QME

rn<-numeric(0)
for(j in 1:nv3){ rn<-c(rn, paste(paste(fac.names[1],':',fac.names[3],sep=''),lf3[j]))}

anavad3<-data.frame("GL"=c(glf3, glE),
"SQ"=c(round(c(SQf3,SQE),5)),
"QM"=c(round(c(QMf3,QME),5)),
"Fc"=c(round(Fcf3,4),''),
"Pr>Fc"=c(round(1-pf(Fcf3,glf3,glE),4),' '))
colnames(anavad3)[5]="Pr>Fc"
rownames(anavad3)=c(rn,"Residuo")
cat('------------------------------------------------------------------------
Quadro da analise de variancia\n------------------------------------------------------------------------\n')
print(anavad3)
cat('------------------------------------------------------------------------\n\n')

ii<-0
for(i in 1:nv3) {
ii<-ii+1
if(1-pf(Fcf3,glf3,glE)[ii]<=sigF){
if(quali[1]==TRUE){
cat('\n\n',fac.names[1],' dentro do nivel ',lf3[i],' de ',fac.names[3],'
------------------------------------------------------------------------')
if(mcomp=='tukey'){
tukey(resp[Fator3==lf3[i]],fatores[,1][Fator3==lf3[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=='duncan'){
duncan(resp[Fator3==lf3[i]],fatores[,1][Fator3==lf3[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='lsd'){
lsd(resp[Fator3==lf3[i]],fatores[,1][Fator3==lf3[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=='lsdb'){
lsdb(resp[Fator3==lf3[i]],fatores[,1][Fator3==lf3[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='sk'){
scottknott(resp[Fator3==lf3[i]],fatores[,1][Fator3==lf3[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='snk'){
snk(resp[Fator3==lf3[i]],fatores[,1][Fator3==lf3[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=="ccboot"){
ccboot(resp[Fator3==lf3[i]],fatores[,1][Fator3==lf3[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=="ccF"){
ccF(resp[Fator3==lf3[i]],fatores[,1][Fator3==lf3[i]],an[9,1],an[9,2],sigT)
 }
}
else{#regressao
cat('\n\n',fac.names[1],' dentro do nivel ',lf3[i],' de ',fac.names[3],'
------------------------------------------------------------------------')
reg.poly(resp[Fator3==lf3[i]], fator1[Fator3==lf3[i]], an[9,1], an[9,2], des3.tab[i+2,1], des3.tab[i+2,2])
}
}
else{cat('\n\n',fac.names[1],' dentro do nivel ',lf3[i],' de ',fac.names[3],'\n')
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp[Fator3==lf3[i]],fatores[,1][Fator3==lf3[i]],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------\n')
}
 }
cat('\n\n')

#Desdobramento de FATOR 3 dentro dos niveis de FATOR 1
cat("\nDesdobrando ", fac.names[3], ' dentro de cada nivel de ', fac.names[1], '
------------------------------------------------------------------------\n')

des4<-aov(resp~Fator1/Fator3)

l3<-vector('list',nv1)
names(l3)<-names(summary(Fator1))
v<-numeric(0)
for(j in 1:nv1) {
for(i in 0:(nv3-2)) v<-cbind(v,i*nv1+j)
l3[[j]]<-v
v<-numeric(0)
}
des4.tab<-summary(des4,split=list('Fator1:Fator3'=l3))[[1]]

#Montando a tabela de ANAVA do des4
glf4=c(as.numeric(des4.tab[3:(nv1+2),1]))

SQf4=c(as.numeric(des4.tab[3:(nv1+2),2]))

QMf4=SQf4/glf4

Fcf4=QMf4/QME

rn<-numeric(0)
for(k in 1:nv1){ rn<-c(rn, paste(paste(fac.names[3],':',fac.names[1],sep=''),lf1[k]))}

anavad4<-data.frame("GL"=c(glf4, glE),
"SQ"=c(round(c(SQf4,SQE),5)),
"QM"=c(round(c(QMf4,QME),5)),
"Fc"=c(round(Fcf4,4),''),
"Pr>Fc"=c(round(1-pf(Fcf4,glf4,glE),4),' '))
colnames(anavad4)[5]="Pr>Fc"
rownames(anavad4)=c(rn,"Residuo")
cat('------------------------------------------------------------------------
Quadro da analise de variancia\n------------------------------------------------------------------------\n')
print(anavad4)
cat('------------------------------------------------------------------------\n\n')

ii<-0
for(i in 1:nv1) {
ii<-ii+1
if(1-pf(Fcf4,glf4,glE)[ii]<=sigF){
if(quali[3]==TRUE){
cat('\n\n',fac.names[3],' dentro do nivel ',lf1[i],' de ',fac.names[1],'
------------------------------------------------------------------------')
if(mcomp=='tukey'){
tukey(resp[Fator1==lf1[i]],fatores[,3][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=='duncan'){
duncan(resp[Fator1==lf1[i]],fatores[,3][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='lsd'){
lsd(resp[Fator1==lf1[i]],fatores[,3][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=='lsdb'){
lsdb(resp[Fator1==lf1[i]],fatores[,3][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='sk'){
scottknott(resp[Fator1==lf1[i]],fatores[,3][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='snk'){
snk(resp[Fator1==lf1[i]],fatores[,3][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=="ccboot"){
ccboot(resp[Fator1==lf1[i]],fatores[,3][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=="ccF"){
ccF(resp[Fator1==lf1[i]],fatores[,3][Fator1==lf1[i]],an[9,1],an[9,2],sigT)
}
}
else{#regressao
cat('\n\n',fac.names[3],' dentro do nivel ',lf1[i],' de ',fac.names[1],'
------------------------------------------------------------------------')
reg.poly(resp[Fator1==lf1[i]], fator3[Fator1==lf1[i]], an[9,1], an[9,2], des4.tab[i+2,1], des4.tab[i+2,2])
}
 }
else{cat('\n\n',fac.names[3],' dentro do nivel ',lf1[i],' de ',fac.names[1],'\n')
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp[Fator1==lf1[i]],fatores[,3][Fator1==lf1[i]],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------\n')
}
}

#Checar o Fator2
if(pvalor[4]>sigF && pvalor[6]>sigF) {
cat('\nAnalisando os efeitos simples do fator ',fac.names[2],'
------------------------------------------------------------------------\n')

i<-2
{
#Para os fatores QUALITATIVOS, teste de Tukey
if(quali[i]==TRUE && pvalor[i]<=sigF) {
cat(fac.names[i])
if(mcomp=='tukey'){
tukey(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
if(mcomp=='duncan'){
duncan(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
if(mcomp=='lsd'){
lsd(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
if(mcomp=='lsdb'){
lsdb(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
if(mcomp=='sk'){
scottknott(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
if(mcomp=='snk'){
snk(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
if(mcomp=="ccboot"){
ccboot(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
if(mcomp=="ccF"){
ccF(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
}

if(quali[i]==TRUE && pvalor[i]>sigF) {
cat(fac.names[i])
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp,fatores[,i],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------')
}

#Para os fatores QUANTITATIVOS, regressao
if(quali[i]==FALSE && pvalor[i]<=sigF){
cat(fac.names[i])
reg.poly(resp, fatores[,i], an[8,1],an[8,2], an[i,1], an[i,2])
}

if(quali[i]==FALSE && pvalor[i]>sigF) {
cat(fac.names[i])
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp,fatores[,i],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------')
}
cat('\n')
}
}
}

#Interacao Fator2*Fator3
if(any(unfold==2.3)){
cat("\n\n\nInteracao",paste(fac.names[2],'*',fac.names[3],sep='')," significativa: desdobrando a interacao
------------------------------------------------------------------------\n")

#Desdobramento de FATOR 2 dentro do niveis de FATOR 3
cat("\nDesdobrando ", fac.names[2], ' dentro de cada nivel de ', fac.names[3], '
------------------------------------------------------------------------\n')

des5<-aov(resp~Fator3/Fator2)

l2<-vector('list',nv3)
names(l2)<-names(summary(Fator3))
v<-numeric(0)
for(j in 1:nv3) {
for(i in 0:(nv2-2)) v<-cbind(v,i*nv3+j)
l2[[j]]<-v
v<-numeric(0)
}
des5.tab<-summary(des5,split=list('Fator3:Fator2'=l2))[[1]]

#Montando a tabela de ANAVA do des5
glf5=c(as.numeric(des5.tab[3:(nv3+2),1]))
SQf5=c(as.numeric(des5.tab[3:(nv3+2),2]))
QMf5=SQf5/glf5
Fcf5=QMf5/QME
rn<-numeric(0)
for(j in 1:nv3){ rn<-c(rn, paste(paste(fac.names[2],':',fac.names[3],sep=''),lf3[j]))}

anavad5<-data.frame("GL"=c(glf5, glE),
"SQ"=c(round(c(SQf5,SQE),5)),
"QM"=c(round(c(QMf5,QME),5)),
"Fc"=c(round(Fcf5,4),''),
"Pr>Fc"=c(round(1-pf(Fcf5,glf5,glE),4),' '))
colnames(anavad5)[5]="Pr>Fc"
rownames(anavad5)=c(rn,"Residuo")
cat('------------------------------------------------------------------------
Quadro da analise de variancia\n------------------------------------------------------------------------\n')
print(anavad5)
cat('------------------------------------------------------------------------\n\n')

ii<-0
for(i in 1:nv3) {
ii<-ii+1
if(1-pf(Fcf5,glf5,glE)[ii]<=sigF){
if(quali[2]==TRUE){
cat('\n\n',fac.names[2],' dentro do nivel ',lf3[i],' de ',fac.names[3],'
------------------------------------------------------------------------')
if(mcomp=='tukey'){tukey(resp[Fator3==lf3[i]],fatores[,2][Fator3==lf3[i]],an[9,1],an[9,2],sigT)}
if(mcomp=='duncan'){duncan(resp[Fator3==lf3[i]],fatores[,2][Fator3==lf3[i]],an[9,1],an[9,2],sigT)}
if(mcomp=='lsd'){lsd(resp[Fator3==lf3[i]],fatores[,2][Fator3==lf3[i]],an[9,1],an[9,2],sigT)}
if(mcomp=='lsdb'){lsdb(resp[Fator3==lf3[i]],fatores[,2][Fator3==lf3[i]],an[9,1],an[9,2],sigT)}
if(mcomp=='sk'){scottknott(resp[Fator3==lf3[i]],fatores[,2][Fator3==lf3[i]],an[9,1],an[9,2],sigT)}
if(mcomp=='snk'){snk(resp[Fator3==lf3[i]],fatores[,2][Fator3==lf3[i]],an[9,1],an[9,2],sigT)}
if(mcomp=="ccboot"){ccboot(resp[Fator3==lf3[i]],fatores[,2][Fator3==lf3[i]],an[9,1],an[9,2],sigT)}
if(mcomp=="ccF"){ccF(resp[Fator3==lf3[i]],fatores[,2][Fator3==lf3[i]],an[9,1],an[9,2],sigT)}
}
else{#regressao
cat('\n\n',fac.names[2],' dentro do nivel ',lf3[i],' de ',fac.names[3],'
------------------------------------------------------------------------')
reg.poly(resp[Fator3==lf3[i]], fator2[Fator3==lf3[i]], an[9,1], an[9,2], des5.tab[i+2,1], des5.tab[i+2,2])
}
}
else{cat('\n\n',fac.names[2],' dentro do nivel ',lf3[i],' de ',fac.names[3],'\n')
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp[Fator3==lf3[i]],fatores[,2][Fator3==lf3[i]],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------\n')
}
}
cat('\n\n')

#Desdobramento de FATOR 3 dentro do niveis de FATOR 2
cat("\nDesdobrando ", fac.names[3], ' dentro de cada nivel de ', fac.names[2], '
------------------------------------------------------------------------\n')

des6<-aov(resp~Fator2/Fator3)

l3<-vector('list',nv2)
names(l3)<-names(summary(Fator2))
v<-numeric(0)
for(j in 1:nv2) {
for(i in 0:(nv3-2)) v<-cbind(v,i*nv2+j)
l3[[j]]<-v
v<-numeric(0)
}
des6.tab<-summary(des6,split=list('Fator2:Fator3'=l3))[[1]]

#Montando a tabela de ANAVA do des6
glf6=c(as.numeric(des6.tab[3:(nv2+2),1]))

SQf6=c(as.numeric(des6.tab[3:(nv2+2),2]))

QMf6=SQf6/glf6

Fcf6=QMf6/QME

rn<-numeric(0)
for(i in 1:nv2){ rn<-c(rn, paste(paste(fac.names[3],':',fac.names[2],sep=''),lf2[i]))}

anavad6<-data.frame("GL"=c(glf6, glE),
"SQ"=c(round(c(SQf6,SQE),5)),
"QM"=c(round(c(QMf6,QME),5)),
"Fc"=c(round(Fcf6,4),''),
"Pr>Fc"=c(round(1-pf(Fcf6,glf6,glE),4),' '))
colnames(anavad6)[5]="Pr>Fc"
rownames(anavad6)=c(rn,"Residuo")
cat('------------------------------------------------------------------------
Quadro da analise de variancia\n------------------------------------------------------------------------\n')
print(anavad6)
cat('------------------------------------------------------------------------\n\n')

ii<-0
for(i in 1:nv2) {
ii<-ii+1
if(1-pf(Fcf6,glf6,glE)[ii]<=sigF){
if(quali[3]==TRUE){
cat('\n\n',fac.names[3],' dentro do nivel ',lf2[i],' de ',fac.names[2],'
------------------------------------------------------------------------')
if(mcomp=='tukey'){
tukey(resp[Fator2==lf2[i]],fatores[,3][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=='duncan'){
duncan(resp[Fator2==lf2[i]],fatores[,3][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='lsd'){
lsd(resp[Fator2==lf2[i]],fatores[,3][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=='lsdb'){
lsdb(resp[Fator2==lf2[i]],fatores[,3][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='sk'){
scottknott(resp[Fator2==lf2[i]],fatores[,3][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='snk'){
snk(resp[Fator2==lf2[i]],fatores[,3][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=="ccboot"){
ccboot(resp[Fator2==lf2[i]],fatores[,3][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=="ccF"){
ccF(resp[Fator2==lf2[i]],fatores[,3][Fator2==lf2[i]],an[9,1],an[9,2],sigT)
 }
}
else{#regressao
cat('\n\n',fac.names[3],' dentro do nivel ',lf2[i],' de ',fac.names[2],'
------------------------------------------------------------------------')
reg.poly(resp[Fator2==lf2[i]], fator3[Fator2==lf2[i]], an[9,1], an[9,2], des6.tab[i+2,1], des6.tab[i+2,2])
}
 }
else{cat('\n\n',fac.names[3],' dentro do nivel ',lf2[i],' de ',fac.names[2],'\n')
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp[Fator2==lf2[i]],fatores[,3][Fator2==lf2[i]],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------\n')
}

}

#Checar o Fator1
if(pvalor[4]>sigF && pvalor[5]>sigF) {
cat('\nAnalisando os efeitos simples do fator ',fac.names[1],'
------------------------------------------------------------------------\n')

i<-1
{
#Para os fatores QUALITATIVOS, teste de Tukey
if(quali[i]==TRUE && pvalor[i]<=sigF) {
cat(fac.names[i])
if(mcomp=='tukey'){
tukey(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
if(mcomp=='duncan'){
duncan(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
if(mcomp=='lsd'){
lsd(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
if(mcomp=='lsdb'){
lsdb(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
if(mcomp=='sk'){
scottknott(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
if(mcomp=='snk'){
snk(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
if(mcomp=="ccboot"){
ccboot(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
if(mcomp=="ccF"){
ccF(resp,fatores[,i],an[8,1],an[8,2],sigT)
}
}

if(quali[i]==TRUE && pvalor[i]>sigF) {
cat(fac.names[i])
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp,fatores[,i],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------')
}

#Para os fatores QUANTITATIVOS, regressao
if(quali[i]==FALSE && pvalor[i]<=sigF){
cat(fac.names[i])
reg.poly(resp, fatores[,i], an[8,1],an[8,2], an[i,1], an[i,2])
}

if(quali[i]==FALSE && pvalor[i]>sigF) {
cat(fac.names[i])
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp,fatores[,i],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------')
}

cat('\n')
}

}
}


#Para interacao tripla significativa, desdobramento
if(any(unfold==3)){
cat("\n\n\nInteracao",paste(fac.names[1],'*',fac.names[2],'*',fac.names[3],sep='')," significativa: desdobrando a interacao
------------------------------------------------------------------------\n")

#Desdobramento de FATOR 1 dentro do niveis de FATOR 2 e do FATOR3
cat("\nDesdobrando ", fac.names[1], ' dentro de cada nivel de ', fac.names[2], 'e',fac.names[3],'
------------------------------------------------------------------------\n')

SQc<-numeric(0)
SQf<-numeric(nv2*nv3)
rn<-numeric(0)

for(i in 1:nv2){
for(j in 1:nv3) {
for(k in 1:nv1) {SQf[(i-1)*nv3+j]=c(SQf[(i-1)*nv3+j]+ sum(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j] & fatores[,1]==lf1[k]])^2) }
rn<-c(rn, paste(paste(fac.names[1],':',sep=''),lf2[i],lf3[j]))
SQc=c(SQc,(sum(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]])^2)/(nv1*J))
}
}
SQf=SQf/J
SQ=SQf-SQc
glf=rep(nv1-1,(nv2*nv3))
QM=SQ/glf
#Montando a tabela de ANAVA do des7

anavad7<-data.frame("GL"=c(glf,glE),
"SQ"=c(SQ,SQE),
"QM"=c(QM,QME),
"Fc"=c(c(round((QM/QME),6)), ' '),
"Pr>Fc"=c(c(round(1-pf(QM/QME,glf,glE),6)),' '))
colnames(anavad7)[5]="Pr>Fc"
rownames(anavad7)=c(rn,"Residuo")
cat('------------------------------------------------------------------------
Quadro da analise de variancia\n------------------------------------------------------------------------\n')
print(anavad7)
cat('------------------------------------------------------------------------\n\n')

ii<-0
for(i in 1:nv2) {
for(j in 1:nv3) {
ii<-ii+1
if(1-pf(QM/QME,glf,glE)[ii]<=sigF){
if(quali[1]==TRUE){
cat('\n\n',fac.names[1],' dentro da combinacao dos niveis ',lf2[i],' de ',fac.names[2],' e ',lf3[j],' de ',fac.names[3],'
------------------------------------------------------------------------')
if(mcomp=='tukey'){tukey(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],an[9,1],an[9,2],sigT)}
if(mcomp=='duncan'){duncan(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],an[9,1],an[9,2],sigT)}
if(mcomp=='lsd'){lsd(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],an[9,1],an[9,2],sigT)}
if(mcomp=='lsdb'){lsdb(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],an[9,1],an[9,2],sigT)}
if(mcomp=='sk'){scottknott(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],an[9,1],an[9,2],sigT)}
if(mcomp=='snk'){snk(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],an[9,1],an[9,2],sigT)}
if(mcomp=="ccboot"){ccboot(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],an[9,1],an[9,2],sigT)}
if(mcomp=="ccF"){ccF(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],an[9,1],an[9,2],sigT)}
}
else{#regressao
cat('\n\n',fac.names[1],' dentro da combinacao dos niveis ',lf2[i],' de ',fac.names[2],' e ',lf3[j],' de ',fac.names[3],'
------------------------------------------------------------------------')
reg.poly(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]], fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]], an[9,1], an[9,2], nv1-1, SQ[ii])
}
}

else{cat('\n\n',fac.names[1],' dentro da combinacao dos niveis ',lf2[i],' de ',fac.names[2],' e ',lf3[j],' de ',fac.names[3],'\n')
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]], fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------\n')
}
}
}
cat('\n\n')

#Desdobramento de FATOR 2 dentro do niveis de FATOR 1 e FATOR 3
cat("\nDesdobrando ", fac.names[2], ' dentro de cada nivel de ', fac.names[1], 'e',fac.names[3],'
------------------------------------------------------------------------\n')

SQc<-numeric(0)
SQf<-numeric(nv1*nv3)
rn<-numeric(0)

for(k in 1:nv1){
for(j in 1:nv3) {
 for(i in 1:nv2) {SQf[(k-1)*nv3+j]=c(SQf[(k-1)*nv3+j]+ sum(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j] & fatores[,2]==lf2[i]])^2) }
rn<-c(rn, paste(paste(fac.names[2],':',sep=''),lf1[k],lf3[j]))
SQc=c(SQc,(sum(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]])^2)/(nv2*J))
}
}
SQf=SQf/J
SQ=SQf-SQc
glf=rep(nv2-1,(nv1*nv3))
QM=SQ/glf
#Montando a tabela de ANAVA do des8


anavad8<-data.frame("GL"=c(glf,glE),
"SQ"=c(SQ,SQE),
"QM"=c(QM,QME),
"Fc"=c(c(round((QM/QME),6)), ' '),
"Pr>Fc"=c(c(round(1-pf(QM/QME,glf,glE),6)),' '))
colnames(anavad8)[5]="Pr>Fc"
rownames(anavad8)=c(rn,"Residuo")
cat('------------------------------------------------------------------------
Quadro da analise de variancia\n------------------------------------------------------------------------\n')
print(anavad8)
cat('------------------------------------------------------------------------\n\n')

ii<-0
for(k in 1:nv1) {
for(j in 1:nv3) {
ii<-ii+1
if(1-pf(QM/QME,glf,glE)[ii]<=sigF){
if(quali[2]==TRUE){
cat('\n\n',fac.names[2],' dentro da combinacao dos niveis ',lf1[k],' de ',fac.names[1],' e ',lf3[j],' de ',fac.names[3],'
------------------------------------------------------------------------')
if(mcomp=='tukey'){
tukey(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],an[9,1],an[9,2],sigT)
}
if(mcomp=='duncan'){
duncan(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='lsd'){
lsd(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],an[9,1],an[9,2],sigT)
}
if(mcomp=='lsdb'){
lsdb(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='sk'){
scottknott(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='snk'){
snk(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],an[9,1],an[9,2],sigT)
}
if(mcomp=="ccboot"){
ccboot(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],an[9,1],an[9,2],sigT)
 }
if(mcomp=="ccF"){
ccF(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],an[9,1],an[9,2],sigT)
 }
}
else{#regressao
cat('\n\n',fac.names[2],' dentro da combinacao dos niveis ',lf1[k],' de ',fac.names[1],' e ',lf3[j],' de ',fac.names[3],'
------------------------------------------------------------------------')
reg.poly(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]], an[9,1], an[9,2], nv2-1, SQ[ii])
}
 }
else{cat('\n\n',fac.names[2],' dentro da combinacao dos niveis ',lf1[k],' de ',fac.names[1],' e ',lf3[j],' de ',fac.names[3],'\n')
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------\n')
}
 }
}

#Desdobramento de FATOR 3 dentro do niveis de FATOR 1 e FATOR 2
cat("\nDesdobrando ", fac.names[3], ' dentro de cada nivel de ', fac.names[1], 'e',fac.names[2],'
------------------------------------------------------------------------\n')

SQc<-numeric(0)
SQf<-numeric(nv1*nv2)
rn<-numeric(0)

for(k in 1:nv1){
for(i in 1:nv2) {
 for(j in 1:nv3) {SQf[(k-1)*nv2+i]=c(SQf[(k-1)*nv2+i]+ sum(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i] & fatores[,3]==lf3[j]])^2) }
rn<-c(rn, paste(paste(fac.names[3],':',sep=''),lf1[k],lf2[i]))
SQc=c(SQc,(sum(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]])^2)/(nv3*J))

}
}
SQf=SQf/J
SQ=SQf-SQc
glf=rep(nv3-1,(nv1*nv2))
QM=SQ/glf
#Montando a tabela de ANAVA do des9


anavad9<-data.frame("GL"=c(glf,glE),
"SQ"=c(SQ,SQE),
"QM"=c(QM,QME),
"Fc"=c(c(round((QM/QME),6)), ' '),
"Pr>Fc"=c(c(round(1-pf(QM/QME,glf,glE),6)),' '))
colnames(anavad9)[5]="Pr>Fc"
rownames(anavad9)=c(rn,"Residuo")
cat('------------------------------------------------------------------------
Quadro da analise de variancia\n------------------------------------------------------------------------\n')
print(anavad9)
cat('------------------------------------------------------------------------\n\n')

ii<-0
for(k in 1:nv1) {
for(i in 1:nv2) {
ii<-ii+1
if(1-pf(QM/QME,glf,glE)[ii]<=sigF){
if(quali[3]==TRUE){
cat('\n\n',fac.names[3],' dentro da combinacao dos niveis ',lf1[k],' de ',fac.names[1],' e ',lf2[i],' de ',fac.names[2],'
------------------------------------------------------------------------')
if(mcomp=='tukey'){
tukey(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=='duncan'){
duncan(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='lsd'){
lsd(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=='lsdb'){
lsdb(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='sk'){
scottknott(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=='snk'){
snk(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],an[9,1],an[9,2],sigT)
}
if(mcomp=="ccboot"){
ccboot(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],an[9,1],an[9,2],sigT)
 }
if(mcomp=="ccF"){
ccF(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],an[9,1],an[9,2],sigT)
 }
}
else{#regressao
cat('\n\n',fac.names[3],' dentro da combinacao dos niveis ',lf1[k],' de ',fac.names[1],' e ',lf2[i],' de ',fac.names[2],'
------------------------------------------------------------------------')
reg.poly(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]], an[9,1], an[9,2], nv3-1, SQ[ii])
}
 }
else{cat('\n\n',fac.names[3],' dentro da combinacao dos niveis ',lf1[k],' de ',fac.names[1],' e ',lf2[i],' de ',fac.names[2],'\n')
cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------\n')
}
 }
}

}
#Saida
out<-list()
out$residuos<-anava$residuals
out$gl.residual<-anava$df.residual
out$coeficientes<-anava$coefficients
out$efeitos<-anava$effects
out$media.Ad<-x
out$valores.ajustados<-anava$fitted.values
out$medias.fator1<-tapply.stat(resp,fatores[,1],mean)
out$medias.fator2<-tapply.stat(resp,fatores[,2],mean)
out$medias.fator3<-tapply.stat(resp,fatores[,3],mean)
tabmedia<-model.tables(aov(resp~Fator1*Fator2*Fator3), "means")
out$medias.dentro12<-tabmedia$tables$`Fator1:Fator2`
out$medias.dentro13<-tabmedia$tables$`Fator1:Fator3`
out$medias.dentro23<-tabmedia$tables$`Fator2:Fator3`
out$medias.dentro123<-tabmedia$tables$`Fator1:Fator2:Fator3`
#if(quali==FALSE && tab[[1]][1,5]<sigF) {out$reg<-reg}
invisible(out)
}
