# ==========================================
#  数据清理模板（Medical Data Cleaning Template）
# ==========================================

# 1. 加载必要包
library(tidyverse)
library(janitor)     # 数据清洗工具，如 clean_names()
library(lubridate)   # 日期处理
library(Hmisc)       # 缺失值工具
library(dplyr)

# 2. 导入数据
data <- read_csv("raw_data.csv")   # 替换为你的文件路径
glimpse(data)                      # 快速查看数据结构

# 3. 统一变量命名
data <- data %>% clean_names()     # 全部变量转为小写、下划线分隔格式

# 4. 检查并处理缺失值 --------------------------------------

# 查看每个变量的缺失比例
na_summary <- sapply(data, function(x) sum(is.na(x)) / length(x))
na_summary <- sort(na_summary, decreasing = TRUE)
print(na_summary)

# 删除缺失过多的变量（例如 >50%）
data <- data[, na_summary < 0.5]

# 用中位数或众数填充
data <- data %>%
  mutate(
    age = ifelse(is.na(age), median(age, na.rm = TRUE), age),
    bmi = ifelse(is.na(bmi), median(bmi, na.rm = TRUE), bmi),
    gender = ifelse(is.na(gender), "Unknown", gender)
  )

# 5. 数据类型转换 -------------------------------------------

# 数值型
num_vars <- c("age", "bmi", "sbp", "dbp")
data[num_vars] <- lapply(data[num_vars], as.numeric)

# 分类变量
cat_vars <- c("gender", "smoking", "diagnosis")
data[cat_vars] <- lapply(data[cat_vars], as.factor)

# 日期变量
data <- data %>%
  mutate(
    admission_date = ymd(admission_date),
    discharge_date = ymd(discharge_date)
  )

# 6. 异常值检测与处理 ---------------------------------------

# 基于IQR检测连续变量的异常值
for (v in num_vars) {
  q <- quantile(data[[v]], probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  lower <- q[1] - 1.5 * iqr
  upper <- q[2] + 1.5 * iqr
  data[[v]][data[[v]] < lower | data[[v]] > upper] <- NA
}

# 检查极端值
summary(data[num_vars])

# 7. 去重与唯一性校验 ---------------------------------------

# 检查重复的患者ID
dup_ids <- data %>% group_by(patient_id) %>% filter(n() > 1)
print(dup_ids)

# 保留第一次记录
data <- data %>% distinct(patient_id, .keep_all = TRUE)

# 8. 逻辑校验 -----------------------------------------------

# 年龄范围合理性
data <- data %>% filter(age >= 0 & age <= 110)

# 检查性别与诊断逻辑（示例：男性不应有“妊娠”类诊断）
if("diagnosis" %in% colnames(data)){
  data <- data %>%
    filter(!(gender == "Male" & str_detect(diagnosis, "Pregnancy")))
}

# 9. 数据合并与特征生成 -------------------------------------

# 示例：计算住院天数
data <- data %>%
  mutate(length_of_stay = as.numeric(discharge_date - admission_date))

# 合并其他表（如患者信息表 + 实验室指标表）
# lab_data <- read_csv("lab_results.csv")
# data <- left_join(data, lab_data, by = "patient_id")

# 10. 导出清洗后的数据 --------------------------------------

write_csv(data, "cleaned_data.csv")
print("✅ 数据清理完成，文件已导出至 cleaned_data.csv")


# 缺失值可视化
library(naniar)
vis_miss(data)  # 直观显示缺失模式
# 异常值分布检查
library(ggplot2)
ggplot(data, aes(x = "", y = bmi)) + geom_boxplot() + theme_minimal()
# 变量分布概览
library(skimr)
skim(data)  # 提供均值、中位数、缺失率、类型统计
