args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
  app_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  app_dir <- normalizePath(".", mustWork = TRUE)
}

gene_set_dir <- file.path(app_dir, "gene_sets")
dir.create(gene_set_dir, showWarnings = FALSE, recursive = TRUE)

obo_url <- "https://purl.obolibrary.org/obo/go/go-basic.obo"
obo_file <- tempfile(fileext = ".obo")
message(sprintf("Downloading GO ontology from %s", obo_url))
utils::download.file(obo_url, obo_file, mode = "wb", quiet = FALSE)

lines <- readLines(obo_file, warn = FALSE, encoding = "UTF-8")

terms <- list()
edges <- list()
current <- NULL

commit_term <- function(term) {
  if (is.null(term) || isTRUE(term$is_obsolete) || is.null(term$id) || !nzchar(term$id) ||
      (!is.null(term$name) && grepl("^obsolete\\b", term$name, ignore.case = TRUE))) {
    return(invisible(NULL))
  }
  if (is.null(term$name) || !nzchar(term$name)) {
    term$name <- term$id
  }
  terms[[length(terms) + 1L]] <<- data.frame(
    term_id = term$id,
    term_name = term$name,
    namespace = if (is.null(term$namespace)) "" else term$namespace,
    stringsAsFactors = FALSE
  )
  if (length(term$edges) > 0) {
    edge_rows <- lapply(term$edges, function(edge) {
      data.frame(
        child_id = term$id,
        parent_id = edge$parent_id,
        relationship = edge$relationship,
        stringsAsFactors = FALSE
      )
    })
    edges[[length(edges) + 1L]] <<- do.call(rbind, edge_rows)
  }
  invisible(NULL)
}

normalize_relation <- function(relation) {
  relation <- trimws(relation)
  relation <- sub("^positively_regulates$", "positive_regulate", relation)
  relation <- sub("^negatively_regulates$", "negative_regulate", relation)
  relation <- sub("^regulates$", "regulate", relation)
  relation
}

for (line in lines) {
  line <- trimws(line)
  if (identical(line, "[Term]")) {
    commit_term(current)
    current <- list(edges = list(), is_obsolete = FALSE)
    next
  }
  if (grepl("^\\[", line)) {
    commit_term(current)
    current <- NULL
    next
  }
  if (is.null(current) || !nzchar(line)) {
    next
  }

  if (startsWith(line, "id: GO:")) {
    current$id <- sub("^id:\\s*", "", line)
  } else if (startsWith(line, "name: ")) {
    current$name <- sub("^name:\\s*", "", line)
  } else if (startsWith(line, "namespace: ")) {
    current$namespace <- sub("^namespace:\\s*", "", line)
  } else if (identical(line, "is_obsolete: true")) {
    current$is_obsolete <- TRUE
  } else if (startsWith(line, "is_a: GO:")) {
    parent_id <- sub("^is_a:\\s*(GO:[0-9]+).*", "\\1", line)
    current$edges[[length(current$edges) + 1L]] <- list(
      parent_id = parent_id,
      relationship = "is_a"
    )
  } else if (startsWith(line, "relationship: ")) {
    relationship <- sub("^relationship:\\s*([^ ]+)\\s+.*", "\\1", line)
    parent_id <- sub("^relationship:\\s*[^ ]+\\s+(GO:[0-9]+).*", "\\1", line)
    relationship <- normalize_relation(relationship)
    if (relationship %in% c("part_of", "regulate", "positive_regulate", "negative_regulate") &&
        grepl("^GO:[0-9]+$", parent_id)) {
      current$edges[[length(current$edges) + 1L]] <- list(
        parent_id = parent_id,
        relationship = relationship
      )
    }
  }
}
commit_term(current)

term_table <- unique(do.call(rbind, terms))
edge_table <- unique(do.call(rbind, edges))
edge_table <- edge_table[
  edge_table$child_id %in% term_table$term_id &
    edge_table$parent_id %in% term_table$term_id,
  ,
  drop = FALSE
]

term_table <- term_table[order(term_table$namespace, term_table$term_id), , drop = FALSE]
edge_table <- edge_table[order(edge_table$child_id, edge_table$relationship, edge_table$parent_id), , drop = FALSE]
rownames(term_table) <- NULL
rownames(edge_table) <- NULL

term_file <- file.path(gene_set_dir, "go_ontology_terms.csv")
edge_file <- file.path(gene_set_dir, "go_ontology_edges.csv")
utils::write.csv(term_table, term_file, row.names = FALSE, na = "")
utils::write.csv(edge_table, edge_file, row.names = FALSE, na = "")

message(sprintf("Wrote %s GO terms to %s", format(nrow(term_table), big.mark = ","), term_file))
message(sprintf("Wrote %s GO edges to %s", format(nrow(edge_table), big.mark = ","), edge_file))
