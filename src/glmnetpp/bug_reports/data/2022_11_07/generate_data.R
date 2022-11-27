load(url("https://github.com/DylanDijk/RepoA/blob/main/reprod_features.rda?raw=true"))
load(url("https://github.com/DylanDijk/RepoA/blob/main/reprod_response.rda?raw=true"))

X = reprod_features
y = reprod_response

write.table(X, "X.csv", sep=',', row.names=F, col.names=F)
write.table(y, "y.csv", sep=',', row.names=F, col.names=F)