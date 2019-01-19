# Install packages (if necessary) and load them.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, 
               caret,
               ModelMetrics,
               scales)

# Set the working directory to the location of the rxnpredict folder.
setwd("C:\\Users\\matt\\Desktop\\rxnpredict\\One-Hot-Encoding_Controls")

# ============================================================================
# Load descriptor and yield data and prepare data for modeling.
# ============================================================================

# Load user-created table containing reaction descriptors.
descriptor.table <- read.csv("R_input-OHE\\descriptor_table-OHE.csv", header=TRUE)

# Scale the descriptor data. Scale parameters are saved in descriptor.data.
descriptor.data <- scale(descriptor.table)
descriptor.scaled <- as.data.frame(descriptor.data)

# Load user-created yield data.
yield.data <- as.numeric(unlist(read.csv("R_input-OHE\\observed_yields-OHE.csv", header=FALSE, stringsAsFactors=FALSE)))

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
saveRDS(rfFit, "R_output-OHE\\rfFit-OHE.rds")

# ============================================================================
# Calculate R^2 and RMSE using test set and generate calibration plot.
# ============================================================================

# Predict yields for test set.
rf.pred <- predict(rfFit, test.scaled)

# Generate *.csv showing predicted and observed yields for the test set (saves to test_set_predicted_yields-OHE.csv).
predicted.yields <- as.data.frame(rf.pred)
predicted.yields$rf.pred <- round(predicted.yields$rf.pred, digits=1)
predicted.yields$yield <- test.scaled$yield
predicted.yields["Error"] <- predicted.yields$yield-predicted.yields$rf.pred
names(predicted.yields)[names(predicted.yields) == 'rf.pred'] <- 'Predicted Yield'
names(predicted.yields)[names(predicted.yields) == 'yield'] <- 'Observed Yield'
write.csv(predicted.yields, 'R_output-OHE\\test_set_predicted_yields-OHE.csv')

# Calculate R^2 and RMS error for test set.
rf.r2 <- cor(rf.pred, test.scaled$yield)^2
rf.rmse <- rmse(rf.pred, test.scaled$yield)

# Calculate R^2 and RMS error for training set (included to compare accuracy of training vs. test sets).
rftrain.pred <- predict(rfFit, training.scaled)
rftrain.r2 <- cor(rftrain.pred, training.scaled$yield)^2
rftrain.rmse <- rmse(rftrain.pred, training.scaled$yield)

# Generate calibration plot of test set (saves to test_set-calibration_plot-OHE.png).
df <- data.frame(x = rf.pred,
                 y = test.scaled$yield)
rsq <- paste(round(rf.r2, digits = 3))
rms <- paste(round(rf.rmse, digits = 1))
p1 <- ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.4) + 
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 101)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(rsq) * "; RMS error = " ~ .(rms) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8)) +
  geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed")
ggsave(file="R_output-OHE\\test_set-calibration_plot-OHE.png", width=5, height=4)

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
ggsave(file="R_output-OHE\\variable_importance_plot-OHE.png", width=8, height=4)

# ============================================================================
# Load descriptors and predict yields for external test set.
# ============================================================================

# Load external test set.
externalset.table <- read.csv("R_input-OHE\\descriptor_table_external_set-OHE.csv", header=TRUE)

# Scale the external set data using the same scaling as for the training and test sets.
externalset.data <- scale(externalset.table,attr(descriptor.data,"scaled:center"),attr(descriptor.data,"scaled:scale"))
externalset.scaled <- as.data.frame(externalset.data)

# Predict yields for external test set.
rf.externalset <- predict(rfFit, externalset.scaled)

# Create table with predicted yields for external test set (saves to external_set_predicted_yields-OHE.csv).
externalset.predictedyields <- as.data.frame(rf.externalset)
externalset.predictedyields $rf.externalset <- round(externalset.predictedyields $rf.externalset, digits=1)
names(externalset.predictedyields )[names(externalset.predictedyields ) == 'rf.externalset'] <- 'Predicted Yield'
write.csv(externalset.predictedyields , 'R_output-OHE\\external_set_predicted_yields-OHE.csv')

# ============================================================================
# Calibration plots for external test substrates with one-hot-encoding.
# ============================================================================

# Load external set observed yields.
external_obs <- as.numeric(unlist(read.csv("R_input-OHE\\external_set_observed_yields-OHE.csv", header=FALSE, stringsAsFactors=FALSE)))

# Store predicted yields for external substrates.
external_pred <- c(rf.externalset)

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

# Generate calibration plot for external substrate 1 (saves to external-calibration_plot-1-OHE.png).
alcohol <- rep(c("a"), times = c(20))
ex1.df <- data.frame(x = external_pred[1:20],
                    y = external_obs[1:20],
                    substrate = alcohol)
ex1rsq <- paste(round(ex1.r2, digits = 3))
ex1rms <- paste(round(ex1.rmse, digits = 1))
ex1.p1 <- ggplot(ex1.df, aes(x = x, y = y, color = "substrate")) +
  geom_point(alpha = 0.6) + 
  scale_color_manual(values=c("black", "red")) +
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(ex1rsq) * "; RMS error = " ~ .(ex1rms) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8), legend.position="none") +
  geom_segment(aes(x=0,xend=100,y=0,yend=100,color="black"), linetype="dashed")
ggsave(file="R_output-OHE\\external-calibration_plot-1-OHE.png", width=5, height=4)

# Generate calibration plot for external substrate 2 (saves to external-calibration_plot-2-OHE.png).
alcohol <- rep(c("a"), times = c(20))
ex2.df <- data.frame(x = external_pred[21:40],
                     y = external_obs[21:40],
                     substrate = alcohol)
ex2rsq <- paste(round(ex2.r2, digits = 3))
ex2rms <- paste(round(ex2.rmse, digits = 1))
ex2.p1 <- ggplot(ex2.df, aes(x = x, y = y, color = "substrate")) +
  geom_point(alpha = 0.6) + 
  scale_color_manual(values=c("black", "darkorange1")) +
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(ex2rsq) * "; RMS error = " ~ .(ex2rms) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8), legend.position="none") +
  geom_segment(aes(x=0,xend=100,y=0,yend=100,color="black"), linetype="dashed")
ggsave(file="R_output-OHE\\external-calibration_plot-2-OHE.png", width=5, height=4)

# Generate calibration plot for external substrate 3 (saves to external-calibration_plot-3-OHE.png).
alcohol <- rep(c("a"), times = c(20))
ex3.df <- data.frame(x = external_pred[41:60],
                     y = external_obs[41:60],
                     substrate = alcohol)
ex3rsq <- paste(round(ex3.r2, digits = 3))
ex3rms <- paste(round(ex3.rmse, digits = 1))
ex3.p1 <- ggplot(ex3.df, aes(x = x, y = y, color = "substrate")) +
  geom_point(alpha = 0.6) + 
  scale_color_manual(values=c("black", "blue")) +
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(ex3rsq) * "; RMS error = " ~ .(ex3rms) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8), legend.position="none") +
  geom_segment(aes(x=0,xend=100,y=0,yend=100,color="black"), linetype="dashed")
ggsave(file="R_output-OHE\\external-calibration_plot-3-OHE.png", width=5, height=4)

# Generate calibration plot for external substrate 4 (saves to external-calibration_plot-4-OHE.png).
alcohol <- rep(c("a"), times = c(20))
ex4.df <- data.frame(x = external_pred[61:80],
                     y = external_obs[61:80],
                     substrate = alcohol)
ex4rsq <- paste(round(ex4.r2, digits = 3))
ex4rms <- paste(round(ex4.rmse, digits = 1))
ex4.p1 <- ggplot(ex4.df, aes(x = x, y = y, color = "substrate")) +
  geom_point(alpha = 0.6) + 
  scale_color_manual(values=c("black", "darkgreen")) +
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(ex4rsq) * "; RMS error = " ~ .(ex4rms) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8), legend.position="none") +
  geom_segment(aes(x=0,xend=100,y=0,yend=100,color="black"), linetype="dashed")
ggsave(file="R_output-OHE\\external-calibration_plot-4-OHE.png", width=5, height=4)

# Generate calibration plot for external substrate 5 (saves to external-calibration_plot-5-OHE.png).
alcohol <- rep(c("a"), times = c(20))
ex5.df <- data.frame(x = external_pred[81:100],
                     y = external_obs[81:100],
                     substrate = alcohol)
ex5rsq <- paste(round(ex5.r2, digits = 3))
ex5rms <- paste(round(ex5.rmse, digits = 1))
ex5.p1 <- ggplot(ex5.df, aes(x = x, y = y, color = "substrate")) +
  geom_point(alpha = 0.6) + 
  scale_color_manual(values=c("black", "purple3")) +
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(ex5rsq) * "; RMS error = " ~ .(ex5rms) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8), legend.position="none") +
  geom_segment(aes(x=0,xend=100,y=0,yend=100,color="black"), linetype="dashed")
ggsave(file="R_output-OHE\\external-calibration_plot-5-OHE.png", width=5, height=4)

# ============================================================================
# Calibration plots for external test substrates with original desciptor model.
# ============================================================================

# Load predicted yields from descriptor model.
external_pred_desc <- c(unlist(read.csv("R_input-OHE\\external_set_predicted_yields-desc.csv", header=FALSE, stringsAsFactors=FALSE)))

# Calculate R^2 and RMSE for external test substrates.
ex1.r2_desc <- cor(external_pred_desc[1:20], external_obs[1:20])^2
ex1.rmse_desc <- rmse(external_pred_desc[1:20], external_obs[1:20])
ex2.r2_desc <- cor(external_pred_desc[21:40], external_obs[21:40])^2
ex2.rmse_desc <- rmse(external_pred_desc[21:40], external_obs[21:40])
ex3.r2_desc <- cor(external_pred_desc[41:60], external_obs[41:60])^2
ex3.rmse_desc <- rmse(external_pred_desc[41:60], external_obs[41:60])
ex4.r2_desc <- cor(external_pred_desc[61:80], external_obs[61:80])^2
ex4.rmse_desc <- rmse(external_pred_desc[61:80], external_obs[61:80])
ex5.r2_desc <- cor(external_pred_desc[81:100], external_obs[81:100])^2
ex5.rmse_desc <- rmse(external_pred_desc[81:100], external_obs[81:100])

# Generate calibration plot for external substrate 1 (saves to external-calibration_plot-1-desc.png).
alcohol <- rep(c("a"), times = c(20))
ex1.df_desc <- data.frame(x = external_pred_desc[1:20],
                     y = external_obs[1:20],
                     substrate = alcohol)
ex1rsq_desc <- paste(round(ex1.r2_desc, digits = 3))
ex1rms_desc <- paste(round(ex1.rmse_desc, digits = 1))
ex1.p1_desc <- ggplot(ex1.df_desc, aes(x = x, y = y, color = "substrate")) +
  geom_point(alpha = 0.6) + 
  scale_color_manual(values=c("black", "red")) +
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(ex1rsq_desc) * "; RMS error = " ~ .(ex1rms_desc) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8), legend.position="none") +
  geom_segment(aes(x=0,xend=100,y=0,yend=100,color="black"), linetype="dashed")
ggsave(file="R_output-OHE\\external-calibration_plot-1-desc.png", width=5, height=4)

# Generate calibration plot for external substrate 2 (saves to external-calibration_plot-2-desc.png).
alcohol <- rep(c("a"), times = c(20))
ex2.df_desc <- data.frame(x = external_pred_desc[21:40],
                     y = external_obs[21:40],
                     substrate = alcohol)
ex2rsq_desc <- paste(round(ex2.r2_desc, digits = 3))
ex2rms_desc <- paste(round(ex2.rmse_desc, digits = 1))
ex2.p1_desc <- ggplot(ex2.df_desc, aes(x = x, y = y, color = "substrate")) +
  geom_point(alpha = 0.6) + 
  scale_color_manual(values=c("black", "darkorange1")) +
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(ex2rsq_desc) * "; RMS error = " ~ .(ex2rms_desc) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8), legend.position="none") +
  geom_segment(aes(x=0,xend=100,y=0,yend=100,color="black"), linetype="dashed")
ggsave(file="R_output-OHE\\external-calibration_plot-2-desc.png", width=5, height=4)

# Generate calibration plot for external substrate 3 (saves to external-calibration_plot-3-desc.png).
alcohol <- rep(c("a"), times = c(20))
ex3.df_desc <- data.frame(x = external_pred_desc[41:60],
                     y = external_obs[41:60],
                     substrate = alcohol)
ex3rsq_desc <- paste(round(ex3.r2_desc, digits = 3))
ex3rms_desc <- paste(round(ex3.rmse_desc, digits = 1))
ex3.p1_desc <- ggplot(ex3.df_desc, aes(x = x, y = y, color = "substrate")) +
  geom_point(alpha = 0.6) + 
  scale_color_manual(values=c("black", "blue")) +
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(ex3rsq_desc) * "; RMS error = " ~ .(ex3rms_desc) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8), legend.position="none") +
  geom_segment(aes(x=0,xend=100,y=0,yend=100,color="black"), linetype="dashed")
ggsave(file="R_output-OHE\\external-calibration_plot-3-desc.png", width=5, height=4)

# Generate calibration plot for external substrate 4 (saves to external-calibration_plot-4-desc.png).
alcohol <- rep(c("a"), times = c(20))
ex4.df_desc <- data.frame(x = external_pred_desc[61:80],
                     y = external_obs[61:80],
                     substrate = alcohol)
ex4rsq_desc <- paste(round(ex4.r2_desc, digits = 3))
ex4rms_desc <- paste(round(ex4.rmse_desc, digits = 1))
ex4.p1_desc <- ggplot(ex4.df_desc, aes(x = x, y = y, color = "substrate")) +
  geom_point(alpha = 0.6) + 
  scale_color_manual(values=c("black", "darkgreen")) +
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(ex4rsq_desc) * "; RMS error = " ~ .(ex4rms_desc) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8), legend.position="none") +
  geom_segment(aes(x=0,xend=100,y=0,yend=100,color="black"), linetype="dashed")
ggsave(file="R_output-OHE\\external-calibration_plot-4-desc.png", width=5, height=4)

# Generate calibration plot for external substrate 5 (saves to external-calibration_plot-5-desc.png).
alcohol <- rep(c("a"), times = c(20))
ex5.df_desc <- data.frame(x = external_pred_desc[81:100],
                     y = external_obs[81:100],
                     substrate = alcohol)
ex5rsq_desc <- paste(round(ex5.r2_desc, digits = 3))
ex5rms_desc <- paste(round(ex5.rmse_desc, digits = 1))
ex5.p1_desc <- ggplot(ex5.df_desc, aes(x = x, y = y, color = "substrate")) +
  geom_point(alpha = 0.6) + 
  scale_color_manual(values=c("black", "purple3")) +
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield', caption = bquote(R^2 ~ " = " ~ .(ex5rsq_desc) * "; RMS error = " ~ .(ex5rms_desc) * "%")) + 
  theme(plot.caption = element_text(hjust = 0.5, size = 8), legend.position="none") +
  geom_segment(aes(x=0,xend=100,y=0,yend=100,color="black"), linetype="dashed")
ggsave(file="R_output-OHE\\external-calibration_plot-5-desc.png", width=5, height=4)

