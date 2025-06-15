library('VIM')
library(pROC)
library(e1071)
library(caret)
library(ggplot2)
library(xgboost)
library(randomForest)
db_detail=readxl::read_excel("D:\\\\database_detail.xlsx")
load(("D:\\totol_dataset.Rda"))
df=metabo_dataset$delta
mutate_data=kNN(df)
mutate_data=t(mutate_data)[1:(ncol(mutate_data)/2),]
idx=sample(nrow(mutate_data),nrow(mutate_data)*0.75)
y=db_detail$Target
y[!(y=="IDH1")]="Not IDH1"
y=y<-as.factor(y)
y_train=y[idx]
x_train=mutate_data[idx,]
y_test=y[-idx]
x_test=mutate_data[-idx,]
BP_result=readxl::read_excel("D:\\BPtest_result.xlsx",col_names = F)$...1
BP_confu=confusionMatrix(factor(BP_result),y_test,positive = 'IDH1')
roc_obj_BP=roc(
  response = y[-idx],
  predictor=readxl::read_excel("D:\\BPvalidation_result.xlsx",col_names = F)$...1[-idx],
  levels = c("Not IDH1", "IDH1"),
  direction = "<"
)
x_matrix <- as.matrix(x_train)
dtrain <- xgb.DMatrix(data = x_matrix, label = as.numeric(y_train) - 1 )
dtest=xgb.DMatrix(data = as.matrix(x_test), label = as.numeric(y_test) - 1 )
params <- list(
  objective = "binary:logistic",  
  eval_metric = "auc",            
  eta = 0.1,                      
  max_depth = 6,                  
  subsample = 0.8,                
  colsample_bytree = 0.8          
)
set.seed(123)
cv_model <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 100,
  nfold = 5,
  stratified = TRUE,
  early_stopping_rounds = 10,
  print_every_n = 10,
  maximize = TRUE
)
best_nrounds <- cv_model$best_iteration
cat("Best iter：", best_nrounds, "\n")
final_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 35,
  watchlist = list(train = dtrain),
  early_stopping_rounds = 10,
  print_every_n = 10,
  maximize = TRUE,
  verbose = 1
)
pred_prob <- predict(final_model, dtest)
pred_label <- ifelse(pred_prob > 0.5, "Not IDH1", "IDH1")
XGboost_confu=confusionMatrix(factor(pred_label),y_test,positive = 'IDH1')
roc_obj_XGBoost <- roc((as.numeric(y_test) - 1), pred_prob,direction='<')
rf_model <- randomForest(
  x = x_train,
  y = y_train,
  ntree = 100,
  mtry = sqrt(ncol(x_train)),
  importance = TRUE,
  proximity = FALSE
)
rf_predic=predict(rf_model, x_test)
rf_confu=confusionMatrix(factor(rf_predic),y_test,positive = 'IDH1')
pred_prob_rf <- predict(rf_model, x_test, type = "prob")[, "IDH1"]
y_train_binary_rf <- ifelse(y_test== "IDH1", 1, 0)
roc_obj_rf <- roc(y_train_binary_rf, pred_prob_rf,direction="<")
svm_model<-svm(x_train, y_train, scale=FALSE, kernel='radial', gamma=5,probability = T)
y_pred <- predict(svm_model, x_test)
probabilities = attr(predict(svm_model, x_test, probability = TRUE), "probabilities")[, "IDH1"]
roc_curve = roc(y_test, probabilities)
plot(roc_curve, col = "darkorange", lwd = 2)
auc_value = auc(roc_curve)
legend("bottomright", legend = paste("AUC =", round(auc_value, 2)), col = "darkorange", lwd = 2)
confusionMatrix(y_pred,y_test,positive = 'IDH1')
roc_curve = roc(y_test, MetTD_predict)
plot(roc_curve, col = "blue", lwd = 2)
auc_value = auc(roc_curve)
legend("bottomright", legend = paste("AUC =", round(auc_value, 2)), col = "blue", lwd = 2)
save(svm_model,file="D:\\svm_model.Rda")
save(idx,file="D:\\train_split_idx.Rda")
MetTD_result=do.call(rbind,lapply(1:111,function(x){
  temp=readxl::read_excel(paste0("D:\\",x,".xlsx"))
  return(data.frame(prediction=ifelse(temp$p[571]<=0.1,"IDH1","Not IDH1"),p=temp$p[571]))
}))
MetTD_predict=do.call(rbind,lapply(1:111,function(x){
  temp=readxl::read_excel(paste0("D:\\",x,".xlsx"))
  return(ifelse(temp$p[571]<=0.1,"IDH1","Not IDH1"))
}))
index=1
roc_obj_MetTD <- roc(response = y[-idx],
               predictor = MetTD_result$p[-idx],
               levels = c("Not IDH1", "IDH1"),
               direction = "<") 
roc_obj_SVM=roc(y_test, attr(predict(svm_model, x_test, probability = TRUE), "probabilities")[, "IDH1"])
plot(roc_obj_XGBoost, main = "ROC - IDH1 classification",
     col = "orange", lwd = 5, cex.lab = 2, font.lab = 7,cex.main = 2)
plot(roc_obj_SVM, col = "red", lwd = 5, add = TRUE)
plot(roc_obj_BP, col = "green", lwd = 5, add = TRUE)
plot(roc_obj_rf, col = "yellow", lwd = 5, add = TRUE)
plot(roc_obj_MetTD, col ="#1c61b6", lwd = 5, add = TRUE)
legend("bottomright",
       legend = c(paste("XGBoost AUC =", round(auc(roc_obj_XGBoost), 5)),
                  paste("SVM AUC =", round(auc(roc_obj_SVM), 5)),
                  paste("BP AUC =", round(auc(roc_obj_BP), 5)),
                  paste("RF AUC =", round(auc(roc_obj_rf), 5)),
                  paste("MetTD AUC =", round(auc(roc_obj_MetTD), 5))
                  ),
       col = c("orange", "red" ,'green','yellow',"blue"
               ),
       cex = 1.5,
       lwd = 3)
ggroc(roc_obj) +
  ggtitle("ROC - IDH1 classification") +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  annotate("text", x = 0.6, y = 0.2, label = paste("AUC =", round(auc_val, 3)),
           color = "#1c61b6", size = 5)
test_indices <- setdiff(1:111, idx)
write.csv(data.frame(train_indices = idx), "D:\\train_indices.csv", row.names = FALSE)
write.csv(data.frame(test_indices = test_indices), "D:\\test_indices.csv", row.names = FALSE)
save.image("D:\\MetTD+SVM+BP+XGBoost+RF.Rda")
