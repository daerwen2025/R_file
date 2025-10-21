# 重新生成“临床试验常用统计分析方法”R脚本模板
# ==========================================
# 临床试验常用统计分析方法总览 (Clinical Trial Statistics Template)
# ==========================================

# 1. 加载常用包 ----------------------------------------------
library(tidyverse)
library(tableone)
library(gtsummary)
library(broom)
library(car)
library(pROC)
library(forestplot)
library(survival)
library(survminer)
library(lme4)
library(DescTools)
library(irr)
library(BlandAltmanLeh)
library(epitools)
library(flextable)

# 2. 导入数据 -------------------------------------------------
# 假设数据文件为 cleaned_data.csv
data <- read_csv("cleaned_data.csv")
glimpse(data)

# ==========================================
# 一、基线比较
# ==========================================

# t检验（连续变量）
t.test(SBP ~ group, data = data)

# 秩和检验（非正态）
wilcox.test(SBP ~ group, data = data)

# 卡方与Fisher检验（分类变量）
chisq.test(table(data$group, data$gender))
fisher.test(table(data$group, data$AE_occurrence))

# ==========================================
# 二、主要终点 - 连续结局 (ANCOVA)
# ==========================================
ancova_model <- lm(Final_BP ~ group + Baseline_BP + age + gender, data = data)
summary(ancova_model)
anova(ancova_model)

# ==========================================
# 三、主要终点 - 二分类结局 (Logistic回归)
# ==========================================
logit_model <- glm(response ~ group + age + gender + BMI, data = data, family = binomial)
summary(logit_model)
exp(cbind(OR = coef(logit_model), confint(logit_model)))

# ==========================================
# 四、生存分析 (Cox + KM)
# ==========================================
fit <- survfit(Surv(time, event) ~ group, data = data)
ggsurvplot(fit, pval = TRUE, risk.table = TRUE)

cox_model <- coxph(Surv(time, event) ~ group + age + gender, data = data)
summary(cox_model)

# ==========================================
# 五、重复测量与混合模型
# ==========================================
# 连续变量的重复测量ANOVA / LMM
lmm <- lmer(SBP ~ time * group + (1|subject_id), data = data)
anova(lmm)
summary(lmm)

# GLMM（Logistic型）
glmm <- glmer(event ~ time * group + (1|subject_id), data = data, family = binomial)
summary(glmm)

# ==========================================
# 六、趋势检验 (剂量-反应关系)
# ==========================================
CochranArmitageTest(table(data$dose_group, data$response))

# ==========================================
# 七、一致性分析
# ==========================================
kappa2(data.frame(rater1 = data$doc1, rater2 = data$doc2))
bland.altman.plot(data$method1, data$method2)

# ==========================================
# 八、亚组与交互分析
# ==========================================
glm_interaction <- glm(response ~ group * gender + age, data = data, family = binomial)
summary(glm_interaction)

# ==========================================
# 九、安全性分析 (AE不良事件)
# ==========================================
table_AE <- table(data$group, data$AE)
chisq.test(table_AE)
riskratio(table_AE)

# ==========================================
# 十、结果导出与报告
# ==========================================
tbl_cox <- tbl_regression(cox_model, exponentiate = TRUE) %>% bold_labels()
as_flex_table(tbl_cox) %>% flextable::save_as_docx(path = "Clinical_Trial_Analysis_Report.docx")

print("✅ 临床试验常用统计分析完成。结果已导出 Clinical_Trial_Analysis_Report.docx。")

