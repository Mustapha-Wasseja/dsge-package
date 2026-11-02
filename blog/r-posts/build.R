# Build the r-posts.com article from dsge-post.Rmd.
#
# Run from this directory:  Rscript build.R
#
# Outputs
#   dsge-post.md          knitted Markdown (code, output and figure links)
#   figures/*.png         the plots, to upload to WordPress
#   wordpress.html        the post body to paste into the WordPress editor
#                         ("Text" tab). R code is wrapped in plain <pre> tags,
#                         as r-posts.com asks, and each figure is replaced by
#                         an [UPLOAD IMAGE: ...] marker showing where to drop
#                         the PNG.
#   preview.html          a standalone page with the figures, for checking
#                         the post before submitting it.
#
# Needs knitr, pandoc and an internet connection (the post downloads FRED
# data).

stopifnot(nzchar(Sys.which("pandoc")))

unlink("figures", recursive = TRUE)
knitr::knit("dsge-post.Rmd", "dsge-post.md", quiet = TRUE)

pandoc_html <- function(md) {
  out <- system2("pandoc",
                 c("-f", "markdown-implicit_figures", "-t", "html",
                   "--no-highlight", "--wrap=none"),
                 input = md, stdout = TRUE)
  paste(out, collapse = "\n")
}

md <- readLines("dsge-post.md", encoding = "UTF-8")

# Drop the YAML header: the title goes into WordPress's title field.
if (length(md) && md[1] == "---") {
  end <- which(md == "---")[2]
  md <- md[-seq_len(end)]
}
html <- pandoc_html(md)

# Plain <pre> blocks for code, as the r-posts.com FAQ asks.
html <- gsub("<pre[^>]*><code[^>]*>", "<pre>", html)
html <- gsub("</code></pre>", "</pre>", html, fixed = TRUE)

# WordPress body: figures become upload markers.
wp <- gsub('<p><img src="figures/([^"]+)"[^>]*/></p>',
           "<p><strong>[UPLOAD IMAGE: \\1]</strong></p>", html)
writeLines(wp, "wordpress.html", useBytes = TRUE)

# Standalone preview page, with the figures embedded so it can be opened
# or shared as a single file.
for (f in list.files("figures", full.names = TRUE)) {
  html <- gsub(paste0('src="', f, '"'),
               paste0('src="', knitr::image_uri(f), '"'), html, fixed = TRUE)
}

title <-"DSGE models in R: estimating a New Keynesian model on US data with dsge"
page <- c(
  "<!doctype html>",
  "<html lang=\"en\"><head><meta charset=\"utf-8\">",
  "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
  paste0("<title>", title, "</title>"),
  "<style>",
  "body{max-width:760px;margin:2rem auto;padding:0 16px;font:17px/1.6 Georgia,serif;color:#222;background:#fff}",
  "h1,h2{font-family:Helvetica,Arial,sans-serif;line-height:1.25}",
  "pre{background:#f5f5f5;padding:12px;overflow-x:auto;font-size:14px;line-height:1.4}",
  "code{font-size:0.9em}",
  "img{max-width:100%;height:auto}",
  "</style></head><body>",
  paste0("<h1>", title, "</h1>"),
  html,
  "</body></html>"
)
writeLines(page, "preview.html", useBytes = TRUE)

message("Wrote dsge-post.md, wordpress.html, preview.html and ",
        length(list.files("figures")), " figures.")
