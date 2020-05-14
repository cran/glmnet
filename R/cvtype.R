cvtype <- function(type.measure="mse",subclass="elnet"){
    type.measures = c("mse","deviance", "class", "auc", "mae","C")
    devname=switch(subclass,
                   elnet="Mean-squared Error",
                   lognet="Binomial Deviance",
                   fishnet="Poisson Deviance",
                   coxnet="Partial Likelihood Deviance",
                   multnet="Multinomial Deviance",
                   mrelnet="Mean-squared Error",
                   glmnetfit="GLM Deviance"
                   )
    typenames = c(deviance = devname, mse = "Mean-Squared Error",
    mae = "Mean Absolute Error",auc = "AUC", class = "Misclassification Error",C="C-index")
    subclass.ch=switch(subclass,
                   elnet=c(1,2,5),
                   lognet=c(2,3,4,1,5),
                   fishnet=c(2,1,5),
                   coxnet=c(2,6),
                   multnet=c(2,3,1,5),
                   mrelnet=c(1,2,5),
                   glmnetfit=c(2,1,5)
                   )
   subclass.type=type.measures[subclass.ch]
   if(type.measure=="default")type.measure=subclass.type[1]
    model.name=switch(subclass,
                   elnet="Gaussian",
                   lognet="Binomial",
                   fishnet="Poisson",
                   coxnet="Cox",
                   multnet="Multinomial",
                   mrelnet="Multi-response Gaussian",
                   glmnetfit="GLM"
                   )
    if(!match(type.measure,subclass.type,FALSE)){
        type.measure=subclass.type[1]
        warning(paste("Only ",paste(subclass.type,collapse=", ")," available as type.measure for ",model.name," models; ", type.measure," used instead",sep=""),call.=FALSE)
    }
    names(type.measure)=typenames[type.measure]
type.measure
    }
