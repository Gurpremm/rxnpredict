# Install packages (if necessary) and load them.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, 
               caret,
               ModelMetrics,
               scales)

# Set the working directory to the location of the rxnpredict folder.
setwd("C:\\Users\\matt\\Desktop\\rxnpredict")

# ============================================================================
# Load descriptor and yield data and prepare data for modeling.
# ============================================================================

# Load user-created table containing reaction descriptors.
descriptor.table <- read.csv("R_input\\descriptor_table.csv", header=TRUE)

# Scale the descriptor data. Scale parameters are saved in descriptor.data.
descriptor.data <- scale(descriptor.table)
descriptor.scaled <- as.data.frame(descriptor.data)

# Load user-created yield data.
yield.data <- as.numeric(unlist(read.csv("R_input\\observed_yields.csv", header=FALSE, stringsAsFactors=FALSE)))

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
train_control <- trainControl(method="cv", number=10, savePredictions=TRUE)

# Train the random forest model.
rfFit <- train(yield ~ ., data=training.scaled, trControl=train_control, method="rf", importance=TRUE)

# Save the trained random forest model.
saveRDS(rfFit, "R_output\\rfFit.rds")

# ============================================================================
# Calculate R^2 and RMSE using test set and generate calibration plot.
# ============================================================================

# Predict yields for test set.
rf.pred <- predict(rfFit, test.scaled)

# Generate *.csv showing predicted and observed yields for the test set (saves to test_set_predicted_yields.csv).
predicted.yields <- as.data.frame(rf.pred)
predicted.yields$rf.pred <- round(predicted.yields$rf.pred, digits=1)
predicted.yields$yield <- test.scaled$yield
predicted.yields["Error"] <- predicted.yields$yield-predicted.yields$rf.pred
names(predicted.yields)[names(predicted.yields) == 'rf.pred'] <- 'Predicted Yield'
names(predicted.yields)[names(predicted.yields) == 'yield'] <- 'Observed Yield'
write.csv(predicted.yields, 'R_output\\test_set_predicted_yields.csv')

# Calculate R^2 and RMS error for test set.
rf.r2 <- cor(rf.pred, test.scaled$yield)^2
rf.rmse <- rmse(rf.pred, test.scaled$yield)

# Calculate R^2 and RMS error for training set (included to compare accuracy of training vs. test sets).
rftrain.pred <- predict(rfFit, training.scaled)
rftrain.r2 <- cor(rftrain.pred, training.scaled$yield)^2
rftrain.rmse <- rmse(rftrain.pred, training.scaled$yield)

# Generate calibration plot of test set (saves to test_set-calibration_plot.png).
df <- data.frame(x = rf.pred,
                 y = test.scaled$yield)
rsq <- paste(round(rf.r2, digits = 3))
rms <- paste(round(rf.rmse, digits = 1))
p1 <- ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.4) + 
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(rsq) * "; RMS error = " ~ .(rms) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8)) +
  geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed")
ggsave(file="R_output\\test_set-calibration_plot.png", width=5, height=4)

# ============================================================================
# Create Variable importance plot.
# ============================================================================

# Read in variable importance from trained rf model.
rf_imp <- importance(rfFit$finalModel)
rf.imp.df <- cbind(as.data.frame(rf_imp), names(rf_imp[, 1]))
colnames(rf.imp.df)[1] <- "IncMSE"
colnames(rf.imp.df)[3] <- "descriptor"

# For descriptor names, replace "_" with " " and "." with "*".
rf.imp.df$descriptor <- gsub("_", " ", rf.imp.df$descriptor)
rf.imp.df$descriptor <- gsub("[.]", "*", rf.imp.df$descriptor)

# Capitalize descriptor names.
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep="", collapse=" ")
}
rf.imp.df$descriptor <- sapply(rf.imp.df$descriptor, simpleCap)

# Plot variable importance (saves to variable_importance_plot.png).
# USER: change '10' on next line to modify minimum percentage cutoff for IncMSE.
p2 <- ggplot(rf.imp.df[rf.imp.df$IncMSE>10, ], aes(x=reorder(descriptor, IncMSE), y=IncMSE)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels = comma) +
  labs(x="", y="Increase in Mean Squared Error (% yield)^2") + 
  coord_flip()
# USER: change 'width' and 'height' parameter on next line to control plot dimensions.
ggsave(file="R_output\\variable_importance_plot.png", width=8, height=4)

# ============================================================================
# Load descriptors and predict yields for external test set.
# ============================================================================

# Load external test set.
externalset.table <- read.csv("R_input\\descriptor_table_external_set.csv", header=TRUE)

# Scale the external set data using the same scaling as for the training and test sets.
externalset.data <- scale(externalset.table,attr(descriptor.data,"scaled:center"),attr(descriptor.data,"scaled:scale"))
externalset.scaled <- as.data.frame(externalset.data)

# Predict yields for external test set.
rf.externalset <- predict(rfFit, externalset.scaled)

# Create table with predicted yields for external test set (saves to external_set_predicted_yields.csv).
externalset.predictedyields <- as.data.frame(rf.externalset)
externalset.predictedyields $rf.externalset <- round(externalset.predictedyields $rf.externalset, digits=1)
names(externalset.predictedyields )[names(externalset.predictedyields ) == 'rf.externalset'] <- 'Predicted Yield'
write.csv(externalset.predictedyields , 'R_output\\external_set_predicted_yields.csv')

# ============================================================================
# Calibration plots for external test substrates.
# ============================================================================

# Load external set observed yields.
external_obs <- as.numeric(unlist(read.csv("R_input\\external_set_observed_yields.csv", header=FALSE, stringsAsFactors=FALSE)))

# Store predicted yields for external substrates.
external_pred <- c(rf.externalset)

# Generate calibration plot for external substrates (saves to external-calibration_plot.png).
ex.r2 <- cor(external_pred, external_obs)^2
ex.rmse <- rmse(external_pred, external_obs)
alcohol <- rep(c("1ag", "1ah", "1ai", "1aj","1ak"), times = c(20,20,20,20,20))
ex.df <- data.frame(x = external_pred,
                    y = external_obs,
                    substrate = alcohol)
exrsq <- paste(round(ex.r2, digits = 3))
exrms <- paste(round(ex.rmse, digits = 1))
ex.p1 <- ggplot(ex.df, aes(x = x, y = y, color = substrate)) +
  scale_color_manual(values=c("red", "darkorange1", "blue", "darkgreen", "purple3", "black")) +
  geom_point(alpha = 0.6) + 
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(exrsq) * "; RMS error = " ~ .(exrms) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8), legend.position="none") +
  geom_segment(aes(x=0,xend=100,y=0,yend=100,color="black"), linetype="dashed")
ggsave(file="R_output\\external-calibration_plot.png", width=5, height=4)

# Calculate R^2 and RMSE for external test substrates.
ex1.r2 <- cor(external_pred[1:20], external_obs[1:20])^2
ex1.rmse <- rmse(external_pred[1:20], external_obs[1:20])
ex2.r2 <- cor(external_pred[21:40], external_obs[21:40])^2
ex2.rmse <- rmse(external_pred[21:40], external_obs[21:40])
ex3.r2 <- cor(external_pred[41:60], external_obs[41:60])^2
ex3.rmse <- rmse(external_pred[41:60], external_obs[41:60])
ex4.r2 <- cor(external_pred[61:80], external_obs[61:80])^2
ex4.rmse <- rmse(external_pred[61:80], external_obs[61:80])
ex5.r2 <- cor(external_pred[81:100], external_obs[81:100])^2
ex5.rmse <- rmse(external_pred[81:100], external_obs[81:100])

# Generate table containing R^2 and RMSE for external test substrates (saves to external_set_stats.csv).
externalsetstats <- matrix(c(1,ex1.r2,ex1.rmse,2,ex2.r2,ex2.rmse,3,ex3.r2,ex3.rmse,4,ex4.r2,ex4.rmse,5,ex5.r2,ex5.rmse),ncol=3,byrow=TRUE)
colnames(externalsetstats) <- c("Ext. Substrate","R^2","RMS Error")
externalsetstats <- as.table(externalsetstats)
write.csv(externalsetstats , 'R_output\\external_set_stats.csv')

# ============================================================================
# Combined calibration plot (used to generate Figure 3 in manuscript).
# ============================================================================

# Generate calibration plot of test set and external substrates (saves to combined-calibration_plot.png).

combined_obs <- c(external_obs,test.scaled$yield)
combined_pred <- c(external_pred,rf.pred)
combined_alcohol <- rep(c("1ag", "1ah", "1ai", "1aj","1ak","test set"), times = c(20,20,20,20,20,192))
combined.df <- data.frame(x = combined_pred,
                          y = combined_obs,
                          substrate = combined_alcohol)
combined.p1 <- ggplot(combined.df, aes(x = x, y = y, color = substrate, shape = substrate, fill = substrate)) +
  scale_color_manual(values=c("#ff3d3d", "#fc9f28", "#415ce2", "#1c9102", "#8841a8", "black", "black")) +
  scale_shape_manual(values=c(22, 24, 21, 25, 23, 20)) +
  scale_fill_manual(values=c("#ff3d3d", "#fc9f28", "#415ce2", "#1c9102", "#8841a8", "black")) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote("Test set:" ~ R^2 ~ " = " ~ .(rsq) * "; RMSE = " ~ .(rms) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8), legend.position="none") +
  geom_segment(aes(x=0,xend=100,y=0,yend=100,color="black"), linetype="dashed")
ggsave(file="R_output\\combined-calibration_plot.png", width=4, height=3.2, dpi=600)