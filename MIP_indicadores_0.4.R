rm(list = ls())

setwd('C:/Users/guilherme.branco/Downloads/MIP - indicadores')

library(xlsx)

# somente um comentário inicial.

# importo a base de dados
mip<-as.data.frame(read.csv2('mip_sudamericana_0.csv',row.names = 1))
legendas<-as.data.frame(read.csv2('legendas.csv',F))

# Defino as Matriz Z, D e X da América do Sul
Z_SA<-mip[1:400,1:400]  
X_SA<-mip[1:400,432]
D_SA<-X_SA - apply(Z_SA,1,sum)

# # Defino as Matriz Z, D e X do Uruguai
# Z_URY<-Z_SA[grepl('URY',names(Z_SA)),grepl('URY',names(Z_SA))]
# X_URY<-X_SA[grepl('URY',names(Z_SA))]
# D_URY<-X_URY - apply(Z_URY,1,sum)



######################################################################
########################## Trabalho 1 ################################
######### Indicador para trás, para frente, upstreamness #############
######################################################################


# Calculo a matriz de Leontief para a América do Sul
A <- as.matrix(Z_SA)%*%solve(diag(X_SA))

# note que é impossível inverter a mátriz diag(X_SA), pois ela é singular.
# Vamos exluir as informações com produção 0. que são:
rownames(Z_SA)[X_SA == 0]

Z_SA_mod<-Z_SA[X_SA!=0,X_SA!=0]
X_SA_mod<-X_SA[X_SA!=0]
D_SA_mod<-X_SA_mod - apply(Z_SA_mod,1,sum)

A_SA <- as.data.frame(as.matrix(Z_SA_mod)%*%solve(diag(X_SA_mod)))
colnames(A_SA) <- colnames(Z_SA_mod)



# defino a matriz inversa de Leontief, L
L_SA <- as.data.frame(solve(diag(1,dim(A_SA))-A_SA))

# Montando os índices de ligação - considerando o próprio setor
L_SA_j = apply(L_SA,2,sum)/length(X_SA_mod)    # somo as linhas de cada coluna
L_SA_i = apply(L_SA,1,sum)/length(X_SA_mod)    # somo as colunas de cada linha
L_SA_ij = sum(L_SA)/(length(X_SA_mod)^2)

# Índices de ligação para trás (via Leontief)
L_tras_comsetor<-L_SA_j/L_SA_ij

# Índices de ligação para frente (via Leontief)
L_frente_comsetor<-L_SA_i/L_SA_ij

# Montando os índices de liLação - sem considerar o próprio setor
L_SA_j_sem = apply(L_SA - diag(diag(as.matrix(L_SA)),dim(L_SA)),2,sum)/(length(X_SA_mod)-1)    # somo as linhas de cada coluna
L_SA_i_sem = apply(L_SA - diag(diag(as.matrix(L_SA)),dim(L_SA)),1,sum)/(length(X_SA_mod)-1)    # somo as colunas de cada linha
L_SA_ij_sem = sum(L_SA - diag(diag(as.matrix(L_SA))))/(length(X_SA_mod)*(length(X_SA_mod)-1) )

# Índices de liLação para trás (via Leontief)
L_tras_semsetor<-L_SA_j_sem/L_SA_ij_sem

# Índices de liLação para frente (via Leontief)
L_frente_semsetor<-L_SA_i_sem/L_SA_ij_sem



# defino a matriz inversa de Gosh, G
G_SA <- as.data.frame(solve(diag(X_SA_mod,dim(L_SA)))%*%as.matrix(L_SA)%*%diag(X_SA_mod,dim(L_SA)))
colnames(G_SA) <- colnames(Z_SA_mod)
rownames(G_SA) <- rownames(Z_SA_mod)

# Montando os índices de ligação - considerando o próprio setor
G_SA_j = apply(G_SA,2,sum)/length(X_SA_mod)    # somo as linhas de cada coluna
G_SA_i = apply(G_SA,1,sum)/length(X_SA_mod)    # somo as colunas de cada linha
G_SA_ij = sum(G_SA)/(length(X_SA_mod)^2)

# Índices de ligação para trás (via Gosh)
G_tras_comsetor<-G_SA_j/G_SA_ij

# Índices de ligação para frente (via Gosh)
G_frente_comsetor<-G_SA_i/G_SA_ij

# Montando os índices de ligação - sem considerar o próprio setor
G_SA_j_sem = apply(G_SA - diag(diag(as.matrix(G_SA)),dim(G_SA)),2,sum)/(length(X_SA_mod)-1)    # somo as linhas de cada coluna
G_SA_i_sem = apply(G_SA - diag(diag(as.matrix(G_SA)),dim(G_SA)),1,sum)/(length(X_SA_mod)-1)    # somo as colunas de cada linha
G_SA_ij_sem = sum(G_SA - diag(diag(as.matrix(G_SA))))/(length(X_SA_mod)*(length(X_SA_mod)-1) )

# Índices de ligação para trás (via Gosh)
G_tras_semsetor<-G_SA_j_sem/G_SA_ij_sem

# Índices de ligação para frente (via Gosh)
G_frente_semsetor<-G_SA_i_sem/G_SA_ij_sem



# Índice de Upstreamness

row.names(mip)

names(mip)

Importpais<-mip[grepl('IM_',row.names(mip)),1:400]
Importag<-apply(Importpais,2,sum)
Importag_mod<-Importag[X_SA!=0]

Exportag<-apply(mip[1:400,grepl('EX_',names(mip))],1,sum)
Exportag_mod<-Exportag[X_SA!=0]


# América do Sul - fechada
# Del<-as.data.frame(t(t(A_SA)/X_SA_mod)*X_SA_mod)

Del2<-(A_SA/rep(X_SA_mod,ncol(A_SA)))*t(rep(X_SA_mod,ncol(A_SA)))


# U_fechada<-solve(diag(1,length(X_SA_mod)) - Del)%*%matrix(1,length(X_SA_mod),1)
U_fechada<-solve(diag(1,nrow(X_SA_mod)) - Del)%*%matrix(1,nrow(X_SA_mod),1)
U_fechada<-solve(diag(1,nrow(X_SA_mod)) - Del2)%*%matrix(1,nrow(X_SA_mod),1)

U_fechada

# América do Sul - aberta
# Del_aberta<-as.data.frame(t(t(A_SA)/(X_SA_mod - Exportag_mod + Importag_mod))*X_SA_mod)

Del_aberta2<-(A_SA/rep(X_SA_mod + Importag_mod - Exportag_mod,ncol(A_SA)))*t(rep(X_SA_mod,ncol(A_SA)))
U_aberta <-solve(diag(1,nrow(X_SA_mod)) - Del_aberta2)%*%matrix(1,nrow(X_SA_mod),1)

U_aberta

# Exportando os resultados

Results<-cbind(L_tras_comsetor,
      L_frente_comsetor,
      L_tras_semsetor,
      L_frente_semsetor,
      G_tras_comsetor,
      G_frente_comsetor,
      G_tras_semsetor,
      G_frente_semsetor,
      U_fechada,
      U_aberta)

Results<-as.data.frame(Results)
names(Results)[9:10]<-c('Up_fechada','Up_aberta')
write.csv2(Results,'Resultados.csv')


# Qual país você deseja?
for (p in PAIS){
unicountry<-Results[grepl(p,row.names(Results)),]
write.csv2(unicountry,paste0('Resultados_',p,'.csv'))
ordem<-apply(unicountry,2,order)
RES<-{}
for (ind in 1:dim(Results)[2]){
  valores<-unicountry[ordem[,ind],ind][c(1:5,(dim(unicountry)[1]-4):dim(unicountry)[1])]
  nomes<-row.names(unicountry)[ordem[,ind]][c(1:5,(dim(unicountry)[1]-4):dim(unicountry)[1])]
  RES[[ind]]<-data.frame(setores = nomes, valores)
}
names(RES)<-names(unicountry)
write.csv2(RES,paste0('indices_',p,'_.csv'))
}



######################################################################
########################## Trabalho 2 ################################
######################## Value Added Share ###########################
######################################################################

# va<-data.frame(VA = as.numeric(mip[row.names(mip)=='VA',]))
# rownames(va)<-colnames(mip)
# va<-subset(va,X_SA != 0 & !is.na(va))

va<-data.frame(1-colSums(A_SA))

V<-va
for (i in 1:9){
V<-cbind(va,V)
}

PAIS<-c('ARG','BRA','BOL','CHL','COL','ECU','PAR','PER','URY','VEN')
for (i in 1:10){
  V[!grepl(row.names(va),pattern = PAIS[i]),i]<-0
}

VAS<-data.frame(t(as.matrix(V))%*%as.matrix(L_SA))
row.names(VAS)<-PAIS
write.xlsx2(VAS,'VAS.xlsx')

# Verifico que a soma das linhas nas colunas é igual a 1.
colSums(VAS)


# Agora vou selecionar dentro dos 3 países:
# os 5 setores com maior valor agregado vindo de dentro do país,
# os 5 com menor valor agregado vindo de dentro, e desses de onde vem a o maior valor agrega.

for (i in 1:10) {
aux<-VAS[grepl(names(VAS),pattern = PAIS[i])][i,]

alto<-aux[order(aux)][(length(aux)-4):length(aux)]       # encontro os maiores VA
baixo<-aux[order(aux)][1:5]                              # menor va
outros<-VAS[names(baixo)][-i,]                           # outros países
va_outro<-sapply(outros,max)                             # maior va que vem do outro país
names(va_outro)<-row.names(outros)[sapply(outros,which.max)]        # nome do país

RES2<-data.frame(Alto_VA = t(alto[1,]),names(alto),Baixo_VA = t(baixo[1,]),names(baixo),va_outro,names(va_outro))
row.names(RES2)<-1:5
names(RES2)<-c('VAS alto','Setor alto','VAS baixo','Setor baixo','VAS baixo vindo do outro pais','outro pais')

write.csv2(RES2,paste0('VAS_',PAIS[i],'_.csv'))
}



# Definindo o vetor E, matriz com as exportações brutas de cada país
Demand_SA<-mip[1:400,401:420]
Demand_SA<-Demand_SA[X_SA!=0,]

E<-matrix(0,nrow = length(Z_SA_mod),ncol = 10)
for (i in 1:10){
# índices para as colunas e linhas
colunas<-!grepl(names(Z_SA_mod),pattern = PAIS[i])
linhas<-grepl(names(Z_SA_mod),pattern = PAIS[i])
# Matriz A e X dados os índices
A_aux<-Z_SA_mod[linhas,colunas]
X_aux<-as.matrix(X_SA_mod[colunas])
# Matriz de demanda dado o índice
colunasd<-!grepl(names(Demand_SA),pattern = PAIS[i])
Demand_aux<-apply(as.matrix(Demand_SA[linhas,colunasd]),1,sum)
E0<-as.matrix(A_aux)%*%X_aux+Demand_aux
E[linhas,i]<-E0
}

VAS_E<-data.frame(t(as.matrix(V))%*%as.matrix(L_SA)%*%E)

VAS_E_per<-matrix(0,nrow=10,ncol = 10)
for (i in 1:10){
VAS_E_per[,i]<-VAS_E[,i]/(apply(VAS_E,2,sum)[i])
}

VAS_E_per<-as.data.frame(VAS_E_per)
colnames(VAS_E_per)<-PAIS
rownames(VAS_E_per)<-PAIS

VAS_E_per

write.csv2(VAS_E_per,'VAS_E_per_.csv')



######################################################################
########################## Trabalho 3 ################################
############################ Grafos ##################################
######################################################################

# carregando pacotes
library(visNetwork)

# Exemplo da apostila da Dai
# criando dados
nodes <- data.frame(id = 1:3, label = c('A','B','C'), color = c('#38B0DE',"#EEAEEE","#66FF66"))
edges <- data.frame(from = c(1,2,3), to = c(2,3,2), arrows = c('to','to','to'))

# esboçar grafo 1 - direcionado
visNetwork(nodes,edges)


####### - Grafo 1 - #######

char1<-{}
for (i in 1:10){
ltraspais<-L_tras_semsetor[grepl(names(L_tras_semsetor),pattern = PAIS[i])]
ind<-which.max(ltraspais)
char1[i]<-names(ltraspais)[ind]
}

char2<-matrix(0,nrow = 3,ncol = 10)
char2<-{}
for (i in 1:10){
zpais<-Z_SA_mod[char1[i]]
char2<-c(char2,rownames(zpais)[order(zpais,decreasing = T)][1:3])
}

char3<-{}
for (i in 1:30){
  zpais<-Z_SA_mod[char2[i]]
  char3<-c(char3,rownames(zpais)[order(zpais,decreasing = T)][1:2])
}

char3

nodes1<-data.frame(id = char1, label = char1, color = '#38B0DE',size = 30,group = 'step 1')
nodes2<-data.frame(id = char2, label = char2, color = '#EEAEEE',size = 25,group = 'step 2')
nodes3<-data.frame(id = char3, label = char3, color = '#66FF66',size = 20,group = 'step 3')

nomes<-c('id','label','color','size','group')
names(nodes1)<-nomes
names(nodes2)<-nomes
names(nodes3)<-nomes

# Seleciono somente os nodes2 que não são iguais aos nodes1
nodes2_sem_rep<-nodes2[!(as.character(nodes2[,1]) %in% as.character(nodes1[,1])),]
# Seleciono somente os nodes3 que não são iguais aos nodes1 e nodes2_sem_rep
nodes3_sem_rep<-nodes3[!(as.character(nodes3[,1]) %in% as.character(rbind(nodes1,nodes2_sem_rep)[,1])),]
# Crio os nós
nodes <- rbind(nodes1,nodes2_sem_rep,nodes3_sem_rep)
# Deleto os nós repetidos
nodes<-nodes[!duplicated(nodes),]

row.names(nodes)<-1:length(nodes[,1])
edges1 <- data.frame(from = rep(char1,each = 3), to = char2, arrows = 'from')
edges2 <- data.frame(from = rep(char2,each = 2), to = char3, arrows = 'from')

edges<-rbind(edges1,edges2)
# Adicionando estilo às setas (pontilhado e cores corretas)
edges<-cbind(edges,color=c(rep('#38B0DE',30),rep('#EEAEEE',60)),dashes = c(rep(F,30),rep(T,60)))

visNetwork(nodes,edges,main = 'Grafos setoriais 1 - América do Sul') %>%
  visGroups(groupname = 'step 1',color = '#38B0DE') %>%
  visGroups(groupname = 'step 2',color = '#EEAEEE') %>%
  visGroups(groupname = 'step 3',color = '#66FF66') %>%
  visLegend()




####### - Grafo 2 - #######

char1<-{}
ind<-order(L_tras_semsetor,decreasing = T)
char1<-names(L_tras_semsetor)[ind][1:10]

char2<-matrix(0,nrow = 3,ncol = 10)
char2<-{}
for (i in 1:10){
  zpais<-Z_SA_mod[char1[i]]
  char2<-c(char2,rownames(zpais)[order(zpais,decreasing = T)][1:3])
}

char3<-{}
for (i in 1:30){
  zpais<-Z_SA_mod[char2[i]]
  char3<-c(char3,rownames(zpais)[order(zpais,decreasing = T)][1:2])
}

nodes1<-data.frame(id = char1, label = char1, color = '#38B0DE',size = 30,group = 'step 1')
nodes2<-data.frame(id = char2, label = char2, color = '#EEAEEE',size = 25,group = 'step 2')
nodes3<-data.frame(id = char3, label = char3, color = '#66FF66',size = 20,group = 'step 3')

nomes<-c('id','label','color','size','group')
names(nodes1)<-nomes
names(nodes2)<-nomes
names(nodes3)<-nomes

# Seleciono somente os nodes2 que não são iguais aos nodes1
nodes2_sem_rep<-nodes2[!(as.character(nodes2[,1]) %in% as.character(nodes1[,1])),]
# Seleciono somente os nodes3 que não são iguais aos nodes1 e nodes2_sem_rep
nodes3_sem_rep<-nodes3[!(as.character(nodes3[,1]) %in% as.character(rbind(nodes1,nodes2_sem_rep)[,1])),]
# Crio os nós
nodes <- rbind(nodes1,nodes2_sem_rep,nodes3_sem_rep)
# Deleto os nós repetidos
nodes<-nodes[!duplicated(nodes),]

row.names(nodes)<-1:length(nodes[,1])
edges1 <- data.frame(from = rep(char1,each = 3), to = char2, arrows = 'from')
edges2 <- data.frame(from = rep(char2,each = 2), to = char3, arrows = 'from')

edges<-rbind(edges1,edges2)
# Adicionando estilo às setas (pontilhado e cores corretas)
edges<-cbind(edges,color=c(rep('#38B0DE',30),rep('#EEAEEE',60)),dashes = c(rep(F,30),rep(T,60)))

visNetwork(nodes,edges,main = 'Grafos setoriais 2 - América do Sul') %>%
  visGroups(groupname = 'step 1',color = '#38B0DE') %>%
  visGroups(groupname = 'step 2',color = '#EEAEEE') %>%
  visGroups(groupname = 'step 3',color = '#66FF66') %>%
  visLegend()


# quais são as pontes existentes nesses grafos?


######################################################################
########################## Trabalho 4 ################################
################## Sofisticação dos serviços  ########################
######################################################################

for (i in 1:10){

# Seleciono quais são os setores de serviço
s36<-paste0(PAIS[i],'s36')
s37<-paste0(PAIS[i],'s37')
s38<-paste0(PAIS[i],'s38')
s39<-paste0(PAIS[i],'s39')
s40<-paste0(PAIS[i],'s40')
# Encontro os índices de cada setor de serviço
ind<-names(D_SA_mod) %in% c(s36,s37,s38,s39,s40)
# Choque de 10% nas demandas de todos os demais setores

## Choque somente nos demais países de determinado país.

D_SA_mod_novo<-D_SA_mod*(1.1)
D_SA_mod_novo[ind]<-D_SA_mod[ind]
# Novo vetor de produção
X_SA_mod_novo<-as.data.frame(as.matrix(L_SA)%*%as.matrix(D_SA_mod_novo))
# Vetores de serviços
X_SA_mod<-as.data.frame(X_SA_mod)
row.names(X_SA_mod)<-row.names(L_SA)

all_s36 <- grepl(names(D_SA_mod),pattern = 's36')
all_s37 <- grepl(names(D_SA_mod),pattern = 's37')
all_s38 <- grepl(names(D_SA_mod),pattern = 's38')
all_s39 <- grepl(names(D_SA_mod),pattern = 's39')
all_s40 <- grepl(names(D_SA_mod),pattern = 's40')
ind1 <- as.logical(all_s36 + all_s37 + all_s38 + all_s39 + all_s40)

ret_serv<-(X_SA_mod_novo[ind1,]-X_SA_mod[ind1,])/X_SA_mod[ind1,]
serv_pais<-ret_serv[i:(i+4)]

# preparo os dados para a regressão
Y<-ret_serv[!(ret_serv %in% serv_pais)]
X<-rep(serv_pais,9)

model<-lm(Y ~ X)
summary(model)
plot(X,Y,main = paste0('Sofisticação setores serviços - ',PAIS[i]))
abline(model,col=2)

}

