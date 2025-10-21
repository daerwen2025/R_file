# ==========================================
# Cox 比例风险回归分析模板
# ==========================================

# 1. 加载包 ----------------------------------------------------
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(forestplot)
library(car)

# 2. 导入数据 --------------------------------------------------
# 假设数据包含以下变量：
# id, followup_time, event (0=生存, 1=死亡), age, gender, BMI, smoking, group
data <- read_csv("cleaned_data.csv")
glimpse(data)

# 3. 创建生存对象 ---------------------------------------------
surv_obj <- Surv(time = data$followup_time, event = data$event)

# 4. 单变量 Cox 回归 ------------------------------------------
univ_vars <- c("age", "gender", "BMI", "smoking", "group")
univ_results <- map_df(univ_vars, function(var){
  f <- as.formula(paste0("surv_obj ~ ", var))
  model <- coxph(f, data = data)
  res <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
  res$variable <- var
  res
})
write_csv(univ_results, "univariate_cox_results.csv")

# 5. 多变量 Cox 回归 ------------------------------------------
multi_model <- coxph(Surv(followup_time, event) ~ age + gender + BMI + smoking + group, data = data)
summary(multi_model)

# 6. 提取 HR 和 95% CI ----------------------------------------
multi_results <- broom::tidy(multi_model, exponentiate = TRUE, conf.int = TRUE)
multi_results <- multi_results %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  rename(HR = estimate, CI_lower = conf.low, CI_upper = conf.high, P = p.value)
write_csv(multi_results, "multivariate_cox_results.csv")

# 7. 检查比例风险假设 -----------------------------------------
cox.zph_test <- cox.zph(multi_model)
print(cox.zph_test)
ggcoxzph(cox.zph_test)  # 可视化检查比例风险假设

# 8. 模型诊断 --------------------------------------------------
# 残差诊断
ggcoxdiagnostics(multi_model, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_minimal())
ggcoxdiagnostics(multi_model, type = "schoenfeld",
                 linear.predictions = FALSE, ggtheme = theme_minimal())

# 9. 绘制森林图 ------------------------------------------------
forest_data <- multi_results %>%
  mutate(term = str_replace_all(term, "_", " ")) %>%
  mutate(HR_label = sprintf("%.2f (%.2f–%.2f)", HR, CI_lower, CI_upper))

tabletext <- cbind(forest_data$term, forest_data$HR_label)
forestplot(labeltext = tabletext,
           mean = forest_data$HR,
           lower = forest_data$CI_lower,
           upper = forest_data$CI_upper,
           xlog = TRUE,
           title = "Hazard Ratios (Cox Model)",
           boxsize = 0.2)

# 10. Kaplan-Meier 生存曲线 -----------------------------------
fit_km <- survfit(Surv(followup_time, event) ~ group, data = data)
ggsurvplot(fit_km,
           data = data,
           pval = TRUE,
           risk.table = TRUE,
           conf.int = TRUE,
           legend.title = "Group",
           palette = c("#2E86C1", "#E74C3C"),
           title = "Kaplan-Meier Survival Curve")

# 11. 输出结果 -------------------------------------------------
write_csv(multi_results, "Cox_Model_Results.csv")
ggsave("KM_Curve.png", width = 6, height = 5, dpi = 300)
print("✅ Cox回归分析完成。结果文件已导出为 Cox_Model_Results.csv 和 KM_Curve.png")


# 检查共线性（VIF）
vif(multi_model)
# 输出表格至 Word
library(gtsummary)
tbl_cox <- tbl_regression(multi_model, exponentiate = TRUE) %>% bold_labels()
as_flex_table(tbl_cox) %>% flextable::save_as_docx(path = "Cox_Model_Table.docx")

# ==========================================
# 医学数据分析 - 线性回归完整流程
# ==========================================

# 1. 加载必要包 ----------------------------------------------
library(tidyverse)
library(broom)          # 模型结果整理
library(car)            # VIF 检查
library(ggplot2)
library(gtsummary)      # 表格输出
library(flextable)

# 2. 导入数据 -------------------------------------------------
# 示例数据包含: SBP（收缩压，连续结局），age, gender, BMI, smoking, group
data <- read_csv("cleaned_data.csv")
glimpse(data)

# 3. 数据预处理 ----------------------------------------------
data <- data %>%
  mutate(
    gender = factor(gender, levels = c("Male", "Female")),
    smoking = factor(smoking, levels = c("No", "Yes")),
    group = factor(group)
  )

# 4. 单变量线性回归 ------------------------------------------
univ_vars <- c("age", "gender", "BMI", "smoking", "group")
univ_results <- map_df(univ_vars, function(var){
  f <- as.formula(paste0("SBP ~ ", var))
  model <- lm(f, data = data)
  res <- broom::tidy(model)
  res$variable <- var
  res
})
write_csv(univ_results, "univariate_linear_results.csv")

# 5. 多变量线性回归 ------------------------------------------
lm_model <- lm(SBP ~ age + gender + BMI + smoking + group, data = data)
summary(lm_model)

# 提取回归系数、置信区间与显著性
lm_results <- broom::tidy(lm_model, conf.int = TRUE)
lm_results <- lm_results %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  rename(Beta = estimate, CI_lower = conf.low, CI_upper = conf.high, P = p.value)
write_csv(lm_results, "multivariate_linear_results.csv")

# 6. 模型诊断 -------------------------------------------------
# 残差正态性
shapiro.test(residuals(lm_model))
# 多重共线性
vif(lm_model)
# 模型诊断图
par(mfrow=c(2,2))
plot(lm_model)

# 7. 预测与残差分析 ------------------------------------------
data$predicted <- predict(lm_model)
data$residuals <- resid(lm_model)

ggplot(data, aes(x = predicted, y = residuals)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, color = "red") +
  theme_minimal() +
  labs(title = "Residual Plot", x = "Predicted SBP", y = "Residuals")

# 8. 模型可视化 ----------------------------------------------
# 散点拟合图（单变量示例）
ggplot(data, aes(x = BMI, y = SBP)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(title = "Relationship between BMI and SBP")

# 回归系数森林图
lm_results <- lm_results %>%
  mutate(term = str_replace_all(term, "_", " "),
         label = sprintf("%.2f (%.2f–%.2f)", Beta, CI_lower, CI_upper))

library(forestplot)
tabletext <- cbind(lm_results$term, lm_results$label)
forestplot(labeltext = tabletext,
           mean = lm_results$Beta,
           lower = lm_results$CI_lower,
           upper = lm_results$CI_upper,
           xlab = "Beta Coefficient (95% CI)",
           title = "Linear Regression Results")

# 9. 模型表格导出（Word） ------------------------------------
tbl_lm <- tbl_regression(lm_model) %>% bold_labels()
as_flex_table(tbl_lm) %>% flextable::save_as_docx(path = "Linear_Regression_Table.docx")

# 10. 结果输出 ------------------------------------------------
write_csv(lm_results, "Linear_Regression_Results.csv")
ggsave("Residuals_Plot.png", width = 6, height = 5, dpi = 300)
print("✅ 线性回归分析完成。结果已输出为 Linear_Regression_Results.csv 与 Linear_Regression_Table.docx。")


# ==========================================
# 医学数据分析 - Logistic回归完整流程
# ==========================================

# 1. 加载必要包 ----------------------------------------------
library(tidyverse)
library(broom)          # 模型结果整理
library(car)            # 多重共线性
library(gtsummary)      # 美观回归表格
library(flextable)
library(pROC)           # ROC曲线
library(forestplot)

# 2. 导入数据 -------------------------------------------------
# 示例：二分类结局 disease (0=无, 1=有)
# 自变量 age, gender, BMI, smoking, group
data <- read_csv("cleaned_data.csv")
glimpse(data)

# 转换变量类型
data <- data %>%
  mutate(
    gender = factor(gender, levels = c("Male", "Female")),
    smoking = factor(smoking, levels = c("No", "Yes")),
    group = factor(group),
    disease = factor(disease, levels = c(0, 1))
  )

# ==========================================
# 3. 单变量 Logistic 回归
# ==========================================
univ_vars <- c("age", "gender", "BMI", "smoking", "group")

univ_results <- map_df(univ_vars, function(var){
  f <- as.formula(paste0("disease ~ ", var))
  model <- glm(f, data = data, family = binomial)
  res <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
  res$variable <- var
  res
})
write_csv(univ_results, "univariate_logistic_results.csv")

# ==========================================
# 4. 多变量 Logistic 回归
# ==========================================
logit_model <- glm(disease ~ age + gender + BMI + smoking + group, data = data, family = binomial)
summary(logit_model)

# OR值与95%CI提取
logit_results <- broom::tidy(logit_model, exponentiate = TRUE, conf.int = TRUE)
logit_results <- logit_results %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  rename(OR = estimate, CI_lower = conf.low, CI_upper = conf.high, P = p.value)
write_csv(logit_results, "multivariate_logistic_results.csv")

# ==========================================
# 5. 模型诊断
# ==========================================

# 多重共线性检查
vif(logit_model)

# Hosmer-Lemeshow 拟合优度检验
library(ResourceSelection)
hoslem.test(data$disease, fitted(logit_model))

# 残差与杠杆值分析
par(mfrow = c(2,2))
plot(logit_model)

# ==========================================
# 6. ROC 曲线与 AUC
# ==========================================
data$pred <- predict(logit_model, type = "response")

roc_obj <- roc(data$disease, data$pred)
plot(roc_obj, print.auc = TRUE, col = "#2E86C1", lwd = 2)
auc(roc_obj)
ggsave("ROC_Curve.png", width = 5, height = 4, dpi = 300)

# ==========================================
# 7. 绘制森林图（OR + 95% CI）
# ==========================================
logit_results <- logit_results %>%
  mutate(term = str_replace_all(term, "_", " "),
         label = sprintf("%.2f (%.2f–%.2f)", OR, CI_lower, CI_upper))

tabletext <- cbind(logit_results$term, logit_results$label)
forestplot(labeltext = tabletext,
           mean = logit_results$OR,
           lower = logit_results$CI_lower,
           upper = logit_results$CI_upper,
           xlog = TRUE,
           title = "Odds Ratios (Logistic Regression)",
           boxsize = 0.2,
           lineheight = "auto")

# ==========================================
# 8. 模型输出（表格与报告）
# ==========================================
tbl_logit <- tbl_regression(logit_model, exponentiate = TRUE) %>% bold_labels()
as_flex_table(tbl_logit) %>% flextable::save_as_docx(path = "Logistic_Regression_Table.docx")

# 导出结果
write_csv(logit_results, "Logistic_Regression_Results.csv")
print("✅ Logistic回归分析完成。结果已输出为 Logistic_Regression_Results.csv 与 Logistic_Regression_Table.docx。")


# 含交互项
interaction_model <- glm(disease ~ age * gender + BMI + smoking, data = data, family = binomial)
summary(interaction_model)

# 分层分析
male_data <- subset(data, gender == "Male")
logit_male <- glm(disease ~ age + BMI + smoking, data = male_data, family = binomial)
summary(logit_male)


# ==========================================
# 医学数据分析 - 广义线性混合模型（GLMM）完整流程
# ==========================================

# 1. 加载必要包 ----------------------------------------------
library(tidyverse)
library(lme4)          # GLMM 核心包
library(lmerTest)      # 提供显著性检验 (p值)
library(broom.mixed)   # 提取模型结果
library(gtsummary)     # 输出美观表格
library(flextable)
library(ggplot2)
library(DHARMa)        # 模型诊断
library(performance)   # 拟合优度与方差成分

# 2. 导入数据 -------------------------------------------------
# 示例数据结构：
# patient_id, hospital_id, disease(0/1), age, gender, BMI, smoking
data <- read_csv("cleaned_data.csv")
glimpse(data)

# 确保因子变量正确
data <- data %>%
  mutate(
    hospital_id = factor(hospital_id),
    gender = factor(gender, levels = c("Male", "Female")),
    smoking = factor(smoking, levels = c("No", "Yes")),
    disease = as.factor(disease)
  )

# ==========================================
# 3. 模型建立（Logistic型 GLMM）
# ==========================================
# 以医院为随机效应（random intercept）
glmm_model <- glmer(disease ~ age + gender + BMI + smoking + (1 | hospital_id),
                    data = data, family = binomial(link = "logit"))
summary(glmm_model)

# 如果为连续结局，可改为 lmer()
# lmixed <- lmer(BMI ~ age + gender + (1|hospital_id), data = data)

# ==========================================
# 4. 提取固定效应与随机效应
# ==========================================
# 固定效应（Fixed effects）
fixef(glmm_model)

# 随机效应（Random effects）
ranef(glmm_model)

# 方差成分
VarCorr(glmm_model)

# ==========================================
# 5. 模型结果整理与导出
# ==========================================
results <- broom.mixed::tidy(glmm_model, conf.int = TRUE, exponentiate = TRUE)
results <- results %>%
  filter(effect == "fixed") %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  rename(OR = estimate, CI_lower = conf.low, CI_upper = conf.high, P = p.value)
write_csv(results, "GLMM_Logistic_Results.csv")

# ==========================================
# 6. 模型诊断
# ==========================================
simulationOutput <- simulateResiduals(fittedModel = glmm_model)
plot(simulationOutput)   # 检查残差正态性与离群点

# 模型性能指标
check_model(glmm_model)
performance::r2_nakagawa(glmm_model)  # R²（marginal + conditional）
icc(glmm_model)                       # 群体间相关系数

# ==========================================
# 7. 可视化
# ==========================================

# 森林图
results <- results %>%
  mutate(term = str_replace_all(term, "_", " "),
         label = sprintf("%.2f (%.2f–%.2f)", OR, CI_lower, CI_upper))

library(forestplot)
tabletext <- cbind(results$term, results$label)
forestplot(labeltext = tabletext,
           mean = results$OR,
           lower = results$CI_lower,
           upper = results$CI_upper,
           xlog = TRUE,
           title = "Odds Ratios (GLMM with Random Intercept)")

# 随机效应可视化（医院层级差异）
library(sjPlot)
plot_model(glmm_model, type = "re", sort.est = TRUE, show.values = TRUE)

# ==========================================
# 8. 模型输出表格（Word）
# ==========================================
tbl_glmm <- tbl_regression(glmm_model, exponentiate = TRUE) %>% bold_labels()
as_flex_table(tbl_glmm) %>% flextable::save_as_docx(path = "GLMM_Regression_Table.docx")

# ==========================================
# 9. 扩展模型：随机斜率模型
# ==========================================
# 如果希望不同医院的年龄效应不同：
glmm_slope <- glmer(disease ~ age + gender + BMI + smoking + (age | hospital_id),
                    data = data, family = binomial)
anova(glmm_model, glmm_slope)  # 模型比较（AIC）

# ==========================================
# 10. 结果总结
# ==========================================
print("✅ GLMM 分析完成。结果文件已导出为 GLMM_Logistic_Results.csv 与 GLMM_Regression_Table.docx。")
