#!/c5/shared/R/2.15.1/bin/Rscript


load("repetitions.RData")
# variable = counts_reduced

# Read all the data
df1= read.table(file="all3_with_first_line_removed_columns.csv",sep=',',header=T)

# the rest of properties
##################################################################
#### Preprocessing of the fingerprints 
##################################################################
d3 = read.table(file="fps_counts.csv",sep=",")



##################################################################
#### Amino Acid Descriptors
#################################################################

##################################################################
#### Amino Acid Descriptors
#################################################################
# Profeat
load("Profeat_descs.RData")
df5=df2
nzv.columns <- nearZeroVar(df5, freqCut = 30/1)
nb_descs = ncol(df5) - length(nzv.columns)
print(nb_descs)
nzv.names <- names(df5)[nzv.columns]
df5 <- df5[, -nzv.columns]
ncol(df5)
#is.highly.correlated <- findCorrelation(cor(df5), 0.95)
#df5 <- df5[,-is.highly.correlated]
ncol(df5)
##
load("zscales_PROTFP_all2.RData")
df2 = seqs_total10
nzv.columns <- nearZeroVar(df2, freqCut = 30/1)
nb_descs = ncol(df2) - length(nzv.columns)
print(nb_descs)
nzv.names <- names(df2)[nzv.columns]
df2 <- df2[, -nzv.columns]
ncol(df2)
#is.highly.correlated <- findCorrelation(cor(df2), 0.95)
#df2 <- df2[,-is.highly.correlated]
ncol(df2)
print("ook")
#################################################################
### Physicochemical
################################################################
load("physchem.RData")
df4 = physchem


print(dim(df1))
print(dim(df2))
print(dim(d3))
print(dim(df4))
print(dim(df5))

###############################################################################################
rows = seq(1,nrow(df1))
df1=cbind(rows,df1)

d3_not_repeated = c()
df2_not_repeated = c()
df4_not_repeated = c()
df5_not_repeated = c()

new_data = c();new_data2 = c();new_data3 = c();new_data4 = c();new_data5 = c()

stdDevs = c()
Bioactivities=c()
ranges=c()
rows_repeated= c()

for (i in 1:nrow(counts_reduced)){
bioactivity=c()
stdDev=c()
datos=c();datos2=c();datos3=c();datos4=c();datos5=c()
subset1=c();subset1_df2=c();subset1_d3=c();subset1_df4=c();subset1_df5=c()
subset2=c();subset2_df2=c();subset2_d3=c();subset2_df4=c();subset2_df5=c()

#first filter
subset1 = df1[which(df1$accession == counts_reduced$V2[i]),]
subset1_df2 = df2[which(df1$accession == counts_reduced$V2[i]),]
subset1_d3 = d3[which(df1$accession == counts_reduced$V2[i]),]
subset1_df4 = df4[which(df1$accession == counts_reduced$V2[i]),]
subset1_df5 = df5[which(df1$accession == counts_reduced$V2[i]),]

#second filter
subset2 = subset1[which(subset1$chembl_id.1 == counts_reduced$V3[i]),]
subset2_df2 = subset1_df2[which(subset1$chembl_id.1 == counts_reduced$V3[i]),]
subset2_d3 = subset1_d3[which(subset1$chembl_id.1 == counts_reduced$V3[i]),]
subset2_df4 = subset1_df4[which(subset1$chembl_id.1 == counts_reduced$V3[i]),]
subset2_df5 = subset1_df5[which(subset1$chembl_id.1 == counts_reduced$V3[i]),]

values = (subset2$standard_value)*10^-9
values = -log(values,base=10)
bioactivity = mean(values)
rangenow=range(values)[2] - range(values)[1]
ranges=c(ranges,rangenow)
#print(range(bioactivity))
stdDev = sd(values)
datos = subset2[1,]
datos2 = subset2_df2[1,]
datos3 = subset2_d3[1,]
datos4 = subset2_df4[1,]
datos5 = subset2_df5[1,]

datos$standard_value = bioactivity
stdDevs = c(stdDevs,stdDev)
Bioactivities=c(Bioactivities,bioactivity)
# We add to the data.frame contatining all the mean values
new_data = rbind(new_data, datos)
new_data2 = rbind(new_data2, datos2)
new_data3 = rbind(new_data3, datos3)
new_data4 = rbind(new_data4, datos4)
new_data5 = rbind(new_data5, datos5)
# We get the rows that which compounds repeated

rows_repeated = c(rows_repeated, subset2$rows)


}

print(length(rows_repeated))

# Data not repeated
df_not_repeated= df1[-rows_repeated,]
d3_not_repeated = d3[-rows_repeated,]
df2_not_repeated = df2[-rows_repeated,]
df4_not_repeated = df4[-rows_repeated,]
df5_not_repeated = df5[-rows_repeated,]

hh= as.vector(df_not_repeated$standard_value)*10^-9
print(range(hh))
hh= -log(hh,base=10)
df_not_repeated$standard_value = hh
print(range(hh))




data_all1 = rbind(df_not_repeated,new_data)
data_all2 = rbind(df2_not_repeated,new_data2)
data_all3 = rbind(d3_not_repeated,new_data3)
data_all4 = rbind(df4_not_repeated,new_data4)
data_all5 = rbind(df5_not_repeated,new_data5)



save(data_all3, file="fps_no_preprocessed.RData")

