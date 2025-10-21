# 创建完整的“医学预测模型合集分析模板”R脚本

# ==========================================
# 医学预测模型合集分析模板 (Predictive Models Template)
# ==========================================

# 1. 加载包 ---------------------------------------------------
library(tidyverse)
library(glmnet)
library(randomForest)
library(e1071)
library(caret)
library(pROC)
library(rms)
library(survival)
library(survminer)
library(xgboost)
library(forestplot)

# 2. 数据导入与预处理 -----------------------------------------
data <- read_csv("cleaned_data.csv")
data <- na.omit(data)
glimpse(data)

x <- model.matrix(event ~ age + gender + BMI + smoking + group, data = data)[, -1]
y <- data$event

# ==========================================
# 一、线性回归模型 --------------------------------------------
lm_model <- lm(outcome ~ age + gender + BMI + smoking, data = data)
summary(lm_model)
data$pred_lm <- predict(lm_model)

# 可视化
ggplot(data, aes(x = pred_lm, y = outcome)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue") +
  theme_minimal() +
  labs(title = "Predicted vs Actual (Linear Model)")

# ==========================================
# 二、Logistic 回归模型 ----------------------------------------
logit_model <- glm(event ~ age + gender + BMI + smoking, data = data, family = binomial)
summary(logit_model)

data$pred_prob <- predict(logit_model, type = "response")
roc_obj <- roc(data$event, data$pred_prob)
plot(roc_obj, print.auc = TRUE, col = "#E74C3C")

# 森林图
logit_results <- broom::tidy(logit_model, exponentiate = TRUE, conf.int = TRUE)
logit_results <- logit_results %>%
  mutate(label = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high))
tabletext <- cbind(logit_results$term, logit_results$label)
forestplot(labeltext = tabletext, mean = logit_results$estimate,
           lower = logit_results$conf.low, upper = logit_results$conf.high,
           xlog = TRUE, title = "Logistic Regression ORs")

# ==========================================
# 三、Cox 比例风险模型 ------------------------------------------
cox_model <- coxph(Surv(time, event) ~ age + gender + BMI, data = data)
summary(cox_model)
ggforest(cox_model, data = data)
fit <- survfit(Surv(time, event) ~ group, data = data)
ggsurvplot(fit, pval = TRUE, risk.table = TRUE)

# ==========================================
# 四、LASSO 回归 (特征选择) -------------------------------------
cvfit <- cv.glmnet(x, y, alpha = 1, family = "binomial")
plot(cvfit)
coef(cvfit, s = "lambda.min")

# ==========================================
# 五、随机森林模型 ----------------------------------------------
rf_model <- randomForest(as.factor(event) ~ age + gender + BMI + smoking + group,
                         data = data, ntree = 500, importance = TRUE)
print(rf_model)
varImpPlot(rf_model, main = "Variable Importance (Random Forest)")

# ==========================================
# 六、支持向量机 (SVM) ------------------------------------------
svm_model <- svm(as.factor(event) ~ age + gender + BMI + smoking, data = data, kernel = "radial")
pred_svm <- predict(svm_model, data)
confusionMatrix(pred_svm, as.factor(data$event))

# ==========================================
# 七、XGBoost 模型 ----------------------------------------------
dtrain <- xgb.DMatrix(data = x, label = y)
params <- list(booster = "gbtree", objective = "binary:logistic",
               eval_metric = "auc", max_depth = 4, eta = 0.1)
xgb_model <- xgb.train(params = params, data = dtrain, nrounds = 200)
xgb.importance(model = xgb_model) %>% xgb.plot.importance()

# ==========================================
# 八、列线图模型 (Nomogram) ------------------------------------
dd <- datadist(data)
options(datadist = 'dd')
nom_model <- lrm(event ~ age + gender + BMI + smoking, data = data, x = TRUE, y = TRUE)
nom <- nomogram(nom_model, fun = plogis, funlabel = "Risk Probability")
plot(nom)

# 校准曲线
cal <- calibrate(nom_model, method = "boot", B = 100)
plot(cal, main = "Calibration Curve")

# ==========================================
# 九、模型验证与性能评估 ---------------------------------------
set.seed(123)
trainIndex <- createDataPartition(data$event, p = .8, list = FALSE)
train <- data[trainIndex, ]
test <- data[-trainIndex, ]

logit_fit <- glm(event ~ age + gender + BMI + smoking, data = train, family = binomial)
test$pred <- predict(logit_fit, newdata = test, type = "response")
roc_test <- roc(test$event, test$pred)
plot(roc_test, col = "#2E86C1", print.auc = TRUE, main = "Model ROC on Test Set")

# ==========================================
# 十、模型导出 --------------------------------------------------
saveRDS(logit_model, "Logistic_Predictive_Model.rds")
saveRDS(xgb_model, "XGBoost_Predictive_Model.rds")

print("✅ 所有预测模型训练、验证与图形输出完成。")

