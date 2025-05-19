# ───────────────────────────────────────────────────────────────────────────────
# 0. Настройка окружения и загрузка пакетов
# ───────────────────────────────────────────────────────────────────────────────
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("cli",  quietly = TRUE)) install.packages("cli")
library(here); library(cli)

init_libs <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(pkgs, update = FALSE, ask = FALSE)
  invisible(lapply(pkgs, library, character.only = TRUE))
}
init_libs(c("edgeR", "ggplot2", "ggrepel"))

# ───────────────────────────────────────────────────────────────────────────────
# 1. Чтение и подготовка данных
# ───────────────────────────────────────────────────────────────────────────────
read_counts <- function(path, sep = "\t") {
  "Читает raw counts, проверяет и округляет при необходимости."
  cli_alert_info("Читаем матрицу счётов из {.file {path}}")
  df <- read.table(here(path),
                   sep        = sep,
                   header     = TRUE,
                   row.names  = 1,
                   check.names= FALSE)
  if (any(df %% 1 != 0, na.rm = TRUE)) {
    cli_alert_warning("Обнаружены дробные значения, округляем…")
    df <- round(df)
  }
  df
}

read_coldata <- function(path) {
  "Читает coldata и проверяет соответствие образцов."
  cli_alert_info("Читаем метаданные из {.file {path}}")
  md <- read.csv(here(path),
                 sep         = ",",
                 header      = TRUE,
                 row.names   = 1,
                 check.names = FALSE)
  md
}

counts  <- read_counts("raw_counts_matrix.csv")
coldata <- read_coldata("coldata.csv")

if (!identical(colnames(counts), rownames(coldata))) {
  stop("Ошибка: имена образцов в counts и coldata не совпадают")
}

# ───────────────────────────────────────────────────────────────────────────────
# 2. Запуск EdgeR-пайплайна
# ───────────────────────────────────────────────────────────────────────────────
run_edgeR <- function(counts,
                      coldata,
                      group_levels,
                      comparisons_list) {
  #' Полный EdgeR-пайплайн с exactTest и сохранением результатов.
  #'
  #' @param counts            матрица raw counts (genes × samples)
  #' @param coldata           data.frame с метаданными (row.names = sample)
  #' @param group_levels      фактор с уровнями в нужном порядке
  #' @param comparisons_list  named list: имя сравнения → вектор из двух уровней c(ref, alt)
  #' @return список таблиц результатов для каждого сравнения
  
  # 1) Создаём DGEList
  group <- factor(coldata$condition, levels = group_levels)
  dge   <- DGEList(counts = counts, group = group, genes = rownames(counts))
  cli_alert_success("Создан DGEList: {nrow(dge)} генов × {ncol(dge)} образцов")
  
  # 2) Фильтрация low-count генов
  keep <- filterByExpr(dge)
  cli_alert_info("Фильтруем: удаляем {sum(!keep)} генов")
  dge   <- dge[keep, , keep.lib.sizes = FALSE]
  
  # 3) TMM-нормализация
  dge <- calcNormFactors(dge)
  cli_alert_success("TMM факторы: {paste(round(dge$samples$norm.factors,3), collapse=', ')}")
  
  # 4) Оценка дисперсии
  dge <- estimateDisp(dge)
  plotBCV(dge, main = "BCV (EdgeR)")
  
  # 5) exactTest и сохранение
  results <- list()
  for (nm in names(comparisons_list)) {
    pair <- comparisons_list[[nm]]
    cli_alert_info("Запускаем exactTest для {nm}")
    
    et  <- exactTest(dge, pair = pair)
    tbl <- topTags(et, n = Inf, sort.by = "PValue")$table
    tbl <- tbl[order(tbl$FDR), ]
    
    results[[nm]] <- tbl
    out_file <- here(paste0("edgeR_", nm, ".csv"))
    write.csv(tbl, file = out_file, row.names = FALSE)
    cli_alert_success("Сохранён результат: {out_file}")
  }
  
  invisible(results)
}

edgeR_results <- run_edgeR(
  counts,
  coldata,
  group_levels = c("Undifferentiated", "Differentiated", "Induced"),
  comparisons  = list(
    Undiff_vs_Diff   = c("Differentiated", "Undifferentiated"),
    Induced_vs_Diff  = c("Differentiated", "Induced")
  )
)

# ───────────────────────────────────────────────────────────────────────────────
# 3. Визуализация результатов EdgeR
# ───────────────────────────────────────────────────────────────────────────────
generate_edgeR_volcano <- function(res, title, outname) {
  "
  Строит и сохраняет вулкано-плот с метками через ggrepel.
  "
  data <- res
  data$Significant <- data$FDR < 0.05 & abs(data$logFC) >= 1
  data$Label       <- ifelse(data$Significant & abs(data$logFC) > 2 & data$FDR < 1e-5,
                             data$genes, "")
  
  p <- ggplot(data, aes(x = logFC, y = -log10(FDR))) +
    geom_point(aes(color = Significant, alpha = Significant), size = 2.5) +
    geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    ggrepel::geom_text_repel(aes(label = Label),
                             size = 4.0,               # крупнее шрифт подписей
                             max.overlaps = 20,
                             box.padding = 0.5,
                             segment.color = "grey50") +
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red3")) +
    scale_alpha_manual(values = c("FALSE" = 0.5,  "TRUE" = 0.8)) +
    labs(title = title, x = "log2 Fold Change", y = "-log10 FDR") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom",
          plot.title       = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major = element_line(size = 0.2))
  
  out_file <- here(paste0(outname, ".png"))
  ggsave(filename = out_file, plot = p, width = 10, height = 7, dpi = 300)
  cli_alert_success("Сохранён вулкано-плот: {out_file}")
}

# создаём и сохраняем оба вулкано-плота
generate_edgeR_volcano(edgeR_results$Undiff_vs_Diff,   "Undifferentiated vs Differentiated", "volcano_Undiff_vs_Diff")
generate_edgeR_volcano(edgeR_results$Induced_vs_Diff,  "Induced vs Differentiated",        "volcano_Induced_vs_Diff")

# ───────────────────────────────────────────────────────────────────────────────
# 4. Финальные сообщения
# ───────────────────────────────────────────────────────────────────────────────

count_edgeR_DEGs <- function(res) {
  sig  <- res$FDR < 0.05
  up   <- sig & res$logFC >= 1
  down <- sig & res$logFC <= -1
  c(total = sum(sig, na.rm = TRUE),
    up    = sum(up,  na.rm = TRUE),
    down  = sum(down, na.rm = TRUE))
}
cli::cli_alert_info("Количество DEG (FDR < 0.05, |logFC| > 1):")
print(count_edgeR_DEGs(edgeR_results$Undiff_vs_Diff))
print(count_edgeR_DEGs(edgeR_results$Induced_vs_Diff))
