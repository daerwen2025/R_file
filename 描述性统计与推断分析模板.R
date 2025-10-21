# ==========================================
# 医学数据分析 - 描述性统计与推断分析模板
# ==========================================

# 加载包 ----------------------------------------------------
library(tidyverse)
library(tableone)
library(gtsummary)
library(compareGroups)
library(FSA)
library(car)
library(janitor)
library(knitr)

# 导入数据 --------------------------------------------------
data <- read_csv("cleaned_data.csv")
glimpse(data)

# ==========================================
# 一、描述性统计
# ==========================================

# 连续变量描述
data %>%
  group_by(group) %>%
  summarise(
    n = n(),
    age_mean = mean(age, na.rm = TRUE),
    age_sd = sd(age, na.rm = TRUE),
    age_median = median(age, na.rm = TRUE),
    age_IQR = IQR(age, na.rm = TRUE)
  )

# 分类变量描述
data %>%
  group_by(group, gender) %>%
  summarise(count = n()) %>%
  mutate(percent = round(count / sum(count) * 100, 1))

# ==========================================
# 二、两组比较
# ==========================================

# t检验（正态）
t.test(age ~ group, data = data, var.equal = TRUE)

# 秩和检验（非参数）
wilcox.test(BMI ~ group, data = data)

# 卡方检验
tbl <- table(data$group, data$smoking)
chisq.test(tbl)

# Fisher精确检验
fisher.test(tbl)

# ==========================================
# 三、多组比较（≥3组）
# ==========================================

# 单因素方差分析 ANOVA
anova_model <- aov(BMI ~ group, data = data)
summary(anova_model)

# 多重比较 Tukey
TukeyHSD(anova_model)

# Kruskal-Wallis 检验
kruskal.test(BMI ~ group, data = data)

# Dunn检验（非参数多重比较）
dunnTest(BMI ~ group, data = data, method = "bonferroni")

# ==========================================
# 四、配对数据分析
# ==========================================

# 配对样本 t检验
t.test(before, after, paired = TRUE)

# Wilcoxon 符号秩检验
wilcox.test(before, after, paired = TRUE)

# McNemar 检验
tbl_paired <- table(before_diagnosis, after_diagnosis)
mcnemar.test(tbl_paired)

# ==========================================
# 五、正态性与方差齐性检验
# ==========================================

shapiro.test(data$age)
bartlett.test(BMI ~ group, data = data)
leveneTest(BMI ~ group, data = data)

# ==========================================
# 六、自动化汇总表格输出
# ==========================================

# compareGroups 自动生成
res <- compareGroups(group ~ age + BMI + gender + smoking, data = data)
createTable(res)
export2word(res, file = "Table1_auto.docx")

# gtsummary + flextable 输出美观表格
tbl <- data %>%
  tbl_summary(by = group,
              statistic = list(all_continuous() ~ "{mean} ± {sd}",
                               all_categorical() ~ "{n} ({p}%)")) %>%
                                                         add_p() %>%
                                                           bold_labels()
                                                         as_flex_table(tbl) %>% flextable::save_as_docx(path = "Descriptive_Comparisons.docx")
  