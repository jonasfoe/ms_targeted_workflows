
library(conflicted)
library(fs)
library(tidyverse)
library(callr)

# root needs to be set to the root of the project for this script to work
# root <- "?"

do_report <- function(type, tag = "", skip_ms1 = FALSE) {
  report_template <- paste0(type, "_Report")
  
  output_temp_base_path <- paste0("!", report_template, "_", tag, "_temp")
  
  qmd_temp_path <- paste0(output_temp_base_path, ".qmd")
  pdf_temp_path <- paste0(output_temp_base_path, ".pdf")
  pdf_final_path <- paste0("!", report_template, "_", tag, ".pdf")
  
  tryCatch(
    {
      file_copy(path(root, "!Reports", paste0("!", report_template, ".qmd")), qmd_temp_path, overwrite = TRUE)
      
      quarto::quarto_render(
        input = qmd_temp_path,
        output_format = "pdf",
        output_file = pdf_temp_path,
        execute = TRUE,
        execute_dir = root,
        execute_params = list(path = getwd(), report_tag = tag, skip_ms1 = skip_ms1),
        cache = FALSE
      )
      
      file_delete(qmd_temp_path)
      message("Placing output at '", pdf_final_path, "'")
      file_copy(pdf_temp_path, pdf_final_path, overwrite = TRUE)
      file_delete(pdf_temp_path)
    },
    finally = {
      c(".qmd", ".pdf", ".aux", ".log", ".tex", ".rmarkdown") %>%
        str_c(qmd_temp_path, .) %>%
        .[file_exists(.)] %>%
        file_delete()
      
      paste0(output_temp_base_path, "_files") %>%
        .[dir_exists(.)] %>%
        dir_delete()
    }
  )
  invisible()
}

list.files() %>%
  str_subset("^!\\w+_Report_.*\\.yaml$") %>%
  as_tibble_col("filename") %>%
  mutate(
    type = str_extract(filename, "(?<=^!)(\\w+)(?=_Report_)"),
    tag = str_extract(filename, "(?<=Report_)(.*)(?=\\.yaml$)")
  ) %>%
  pwalk(function(type, tag, ...) do_report(type, tag))
