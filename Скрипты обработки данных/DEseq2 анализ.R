# ───────────────────────────────────────────────────────────────────────────────
# 0. Настройка окружения и загрузка пакетов
# ───────────────────────────────────────────────────────────────────────────────
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("cli",  quietly = TRUE)) install.packages("cli")
library(here)
library(cli)

# Функция инициализации Bioconductor-пакетов
init_libs <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(pkgs, update = FALSE, ask = FALSE)
  invisible(lapply(pkgs, library, character.only = TRUE))
}

# Устанавливем и загружаем нужные пакеты
init_libs(c("DESeq2", "apeglm", "EnhancedVolcano", "ggplot2", "cli"))

# ───────────────────────────────────────────────────────────────────────────────
# 1. Чтение и подготовка данных
# ───────────────────────────────────────────────────────────────────────────────
read_counts <- function(path, sep = "\t") {
  cli::cli_alert_info("Читаем матрицу raw counts из {.file {path}}")
  df <- read.table(here(path),
                   sep       = sep,
                   header    = TRUE,
                   row.names = 1,
                   check.names = FALSE)
  if (!all(round(df) == df, na.rm = TRUE)) {
    cli::cli_alert_warning("Обнаружены дробные счёты, выполняем округление")
    df <- round(df)
  }
  df
}

read_coldata <- function(path) {
  cli::cli_alert_info("Читаем coldata из {.file {path}}")
  df <- read.csv(here(path),
                 sep        = ",",
                 header     = TRUE,
                 row.names  = 1,
                 check.names = FALSE)
  df
}

counts  <- read_counts("raw_counts_matrix.csv")
coldata <- read_coldata("coldata.csv")

# Проверяем соответствие образцов
if (!identical(rownames(coldata), colnames(counts))) {
  stop("Error: Образцы в counts и coldata не совпадают!")
}

# ───────────────────────────────────────────────────────────────────────────────
# 2. Подготовка DESeq2-объекта
# ───────────────────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~condition
)

# Фильтрация: минимум 10 ридов в ≥3 образцах
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cli::cli_alert_success("После фильтрации осталось {nrow(dds)} генов")

# Устанавливаем базовый уровень
dds$condition <- relevel(dds$condition, ref = "Differentiated")

# ───────────────────────────────────────────────────────────────────────────────
# 3. Запуск DESeq2 и извлечение результатов с LFC-shrinkage
# ───────────────────────────────────────────────────────────────────────────────
dds <- DESeq(dds)
resultsNames(dds)

get_shrunk <- function(coef_name, suffix) {
  cli::cli_alert_info("Выполняем LFC-shrinkage для {coef_name}")
  res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  out <- as.data.frame(res[order(res$padj), ])
  file <- here(paste0("DESeq2_", suffix, ".csv"))
  write.csv(out, file, na = "", row.names = TRUE)
  cli::cli_alert_success("Сохранён {file}")
  res
}

res_undiff_vs_diff   <- get_shrunk("condition_Undifferentiated_vs_Differentiated",  "Undiff_vs_Diff")
res_induced_vs_diff  <- get_shrunk("condition_Induced_vs_Differentiated",     "Induced_vs_Diff")

# ───────────────────────────────────────────────────────────────────────────────
# 4. Визуализация
# ───────────────────────────────────────────────────────────────────────────────
# 4.1 PCA
vsd <- vst(dds, blind = TRUE)
pca_plot <- plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA: условные группы") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12)
  )
ggsave(here("PCA_plot.png"), pca_plot, width = 8, height = 6, dpi = 300)
cli::cli_alert_success("Сохранён PCA_plot.png")

# 4.2 Функция для вулкано-плота с регулировкой шрифтов
make_volcano <- function(res, subtitle, outname) {
  # Параметры по умолчанию
  p <- EnhancedVolcano(
    res,
    lab            = rownames(res),
    x              = "log2FoldChange",
    y              = "padj",
    pCutoff        = 0.05,
    FCcutoff       = 1,
    title          = NULL,
    subtitle       = subtitle,
    caption        = paste("Total points:", nrow(res)),
    legendPosition = "bottom",
    pointSize      = 3.5,
    labSize        = 5.0,       # увеличенный шрифт подписей точек
    colAlpha       = 0.8,
    shape          = 19,
    max.overlaps   = 30,
    boxedLabels    = TRUE,
    drawConnectors = TRUE,
    widthConnectors= 0.5,
    col            = c("grey30", "forestgreen", "royalblue", "red3"),
    gridlines.major= FALSE,
    gridlines.minor= FALSE,
    border         = "partial",
    borderWidth    = 1.5,
    borderColour   = "black"
  ) +
    theme_classic() +
    theme(
      plot.title  = element_blank(),
      plot.margin = margin(t = 40, unit = "pt"),
      legend.text = element_text(size = 12),        # чуть больше шрифт легенды
      axis.title  = element_text(size = 14),
      axis.text   = element_text(size = 12)
    ) +
    # Аннотации в углах
    annotate(
      "text",
      x     = -Inf, y     =  Inf,
      label = "Дифференцированные",
      hjust = -0.1, vjust = 2,
      size  = 5,
      color = "black"
    ) +
    annotate(
      "text",
      x     =  Inf, y     =  Inf,
      label = "Недифференцированные",
      hjust =  1.1, vjust = 2,
      size  = 5,
      color = "black"
    )
  
  # Сохраняем
  ggsave(
    filename = here::here(paste0(outname, ".png")),
    plot     = p,
    width    = 10,
    height   = 8,
    dpi      = 300
  )
  cli::cli_alert_success("Сохранён {outname}.png")
}

make_volcano(res_undiff_vs_diff,  "Undifferentiated vs Differentiated", "volcano_Undiff_vs_Diff")
make_volcano(res_induced_vs_diff, "Induced vs Differentiated",       "volcano_Induced_vs_Diff")

# ───────────────────────────────────────────────────────────────────────────────
# 5. Финальные сообщения
# ───────────────────────────────────────────────────────────────────────────────
cli::cli_alert_info("Количество DEG (padj<0.05, |log2FC|>1):")
count_DEGs <- function(res) {
  sig   <- res$padj < 0.05
  up    <- sig & (res$log2FoldChange >=  1)
  down  <- sig & (res$log2FoldChange <= -1)
  c(total = sum(sig, na.rm = TRUE),
    up    = sum(up,  na.rm = TRUE),
    down  = sum(down,na.rm = TRUE))
}
print(count_DEGs(res_undiff_vs_diff))
print(count_DEGs(res_induced_vs_diff))