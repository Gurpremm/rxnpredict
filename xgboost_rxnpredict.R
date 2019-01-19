
library(tidyverse)
library(xgboost)
pacman::p_load(ggplot2, 
               ModelMetrics,
               scales)
pacman::p_load(caret,randomForest)

# Set the working directory to the location of the rxnpredict folder.
setwd("/Users/gurprem/Downloads/rxnpredict")

# ============================================================================
# Load descriptor and yield data and prepare data for modeling.
# ============================================================================

# Load user-created table containing reaction descriptors.
descriptor.table <- read.csv("R_input/descriptor_table.csv", header=TRUE)

# Scale the descriptor data. Scale parameters are saved in descriptor.data.
descriptor.data <- scale(descriptor.table)
descriptor.scaled <- as.data.frame(descriptor.data)

# Load user-created yield data.
yield.data <- as.numeric(unlist(read.csv("R_input/observed_yields.csv", header=FALSE, stringsAsFactors=FALSE)))

# Append the yield data to the descriptor table.
descriptor.scaled$yield <- yield.data

# ============================================================================
# Split data and train random forest model.
# ============================================================================

# Split into training and test set (70/30).
set.seed(1751)
size <- round(0.70*nrow(descriptor.scaled))
training <- sample(nrow(descriptor.scaled), size=size, replace=FALSE)
training.scaled <- descriptor.scaled[training,]
test.scaled <- descriptor.scaled[-training,]

# 10-fold cross-validation.
#train_control <- trainControl(method="cv", number=10, savePredictions=TRUE)
train_contro <- trainControl(method="repeatedcv", repeats = 1,number = 3)

#dmatrix form
#labels
train_y = training.scaled[,'yield']
train_x = training.scaled[, names(training.scaled) !='yield']

test_y = test.scaled[,'yield']
test_x = test.scaled[, names(test.scaled) !='yield']

#dmatrix for test and train
dtrain = xgb.DMatrix(data =  as.matrix(train_x),label = train_y)
dtest = xgb.DMatrix(data =  as.matrix(test_x),label = test_y)

watchlist = list(train=dtrain, test=dtest)

# Train the gradient boost model.
model <- xgb.train(data=dtrain,max.depth=5,eta = 0.3, gamma=5,nround = 110,
                   watchlist = watchlist,objective = "reg:linear",early_stopping_rounds = 60,lambda=0,alpha=1,
                   min_child_weight= 2,print_every_n = 10,trControl = train_control,verbose = 1)

#model$bestTune
#xgb.save(model, "/Users/gurprem/Desktop/xgb.model")
#xgb.plot.tree(model = model)

# Make predictions on the test data
train_predictions <-model %>% predict(dtrain)
predictions <- model %>% predict(dtest)
#head(predictions)
#print(predictions)

#write the predcited outcomes to file
#write.csv(train_predictions, file = "/Users/gurprem/Desktop/train_predictions.csv")
#write.csv(predictions,file="/Users/gurprem/Desktop/test_predictions.csv")

#write.csv(train_y,file="/Users/gurprem/Desktop/train_actual_yield.csv")
#write.csv(test_y,file="/Users/gurprem/Desktop/test_actual_yield.csv")

# Compute the average prediction error RMSE
train_rmse <- rmse(train_predictions, training.scaled$yield)
train_r2 <- cor(train_predictions, training.scaled$yield)^2
test_rmse <- RMSE(predictions, test.scaled$yield)
test_r2 <- cor(predictions, test.scaled$yield)^2


# Load external test set.
externalset.table <- read.csv("R_input/descriptor_table_external_set.csv", header=TRUE)

# Scale the external set data using the same scaling as for the training and test sets.
externalset.data <- scale(externalset.table,attr(descriptor.data,"scaled:center"),attr(descriptor.data,"scaled:scale"))
externalset.scaled <- as.data.frame(externalset.data)

# Load external set observed yields.
external_obs <- as.numeric(unlist(read.csv("R_input/external_set_observed_yields.csv", header=FALSE, stringsAsFactors=FALSE)))
external_obs.data <- scale(externalset.table,attr(descriptor.data,"scaled:center"),attr(descriptor.data,"scaled:scale"))
external_obs.scaled<- as.data.frame(external_obs.data)

# Append the yield data to the external descriptor table.
externalset.scaled$yield <- external_obs

#dbmatrix form external set
ext_y = externalset.scaled[,'yield']
ext_x = externalset.scaled[, names(externalset.scaled) !='yield']

#etest_y = test.scaled[,'yield']
#etest_x = test.scaled[, names(test.scaled) !='yield']
ed_ext = xgb.DMatrix(data =  as.matrix(ext_x),label = ext_y)
#edtest = xgb.DMatrix(data =  as.matrix(test.scaled),label = ext_y)
# Predict yields for external test set.
ext_predictions <- model %>% predict(ed_ext)

# Store predicted yields for external set 
#write.csv(ext_predictions,file="/Users/gurprem/Desktop/external_predictions.csv")

# Generate calibration plot for external substrates (saves to external-calibration_plot.png).
ex.r2 <- cor(ext_predictions, external_obs)^2
ex.rmse <- RMSE(ext_predictions, external_obs)

exgb.df <- data.frame(x = ext_predictions,
                      y = external_obs)
varImp(model)
# Generate calibration plot of test set (saves to test_set-calibration_plot.png).
df_test_plot <- data.frame(x = predictions,
                           y = test.scaled$yield)
#rsq <- paste(round(rf.r2, digits = 3))
#rms <- paste(round(rf.rmse, digits = 1))
p1 <- ggplot(df_test_plot, aes(x = x, y = y)) +
  geom_point(color='steelblue',alpha = 0.4) + 
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(test_r2) * "; RMS error = " ~ .(test_rmse) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8)) +
  geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed")
ggsave(file="/Users/gurprem/Desktop/test_set_plot.png", width=5, height=4)

#plot for external set
ex.p1 <- ggplot(exgb.df, aes(x = x, y = y)) +
  scale_color_manual(values=c("red", "darkorange1", "blue", "darkgreen", "purple3", "black")) +
  geom_point(alpha = 0.6) + 
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(ex.r2) * "; RMS error = " ~ .(ex.rmse) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8), legend.position="none") +
  geom_segment(aes(x=0,xend=100,y=0,yend=100,color="black"), linetype="dashed")
ggsave(file="R_output\\external-calibration_plot_test_001.png", width=5, height=4)

mat <- xgb.importance(feature_names = colnames(descriptor.data),model = model)
xgb.plot.importance(importance_matrix = mat[1:20]) #first 20 variables
print(mat)

summary(model)


rfFit <- train(yield ~ ., data=training.scaled,trControl=train_control, method="rf", importance=TRUE)
rf.pred <- predict(rfFit, test.scaled)
rf.r2 <- cor(rf.pred, test.scaled$yield)^2
rf.rmse <- rmse(rf.pred, test.scaled$yield)
rf.external <- predict(rfFit, externalset.scaled)
rfex.r2 <- cor(rf.external, external_obs)^2
rfex.rmse <- rmse(rf.external, external_obs)


