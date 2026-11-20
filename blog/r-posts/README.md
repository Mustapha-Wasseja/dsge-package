# Guest post for r-posts.com / R-bloggers

Source and build files for a post introducing dsge, to be submitted at
<https://r-posts.com/>. Approved posts appear on R-bloggers.

| File | What it is |
|------|------------|
| `dsge-post.Rmd` | The post. Edit this. |
| `build.R` | Knits the post and writes the files below (`Rscript build.R`, about a minute; downloads FRED data). |
| `wordpress.html` | Post body to paste into WordPress. |
| `figures/*.png` | The plots to upload. |
| `preview.html` | Standalone preview with the plots embedded. |
| `dsge-post.md` | Knitted Markdown, for review on GitHub. |

The post uses features new in dsge 1.2.1, including `autoplot()`, the
`response` argument of the IRF plot, and the corrected `irf_2nd_order()`.
Submit it once 1.2.1 is on CRAN, so `install.packages("dsge")` gives readers
a version that runs every snippet. Rebuild before submitting so the numbers
come from the released version.

## Submitting

r-posts.com rules (from its FAQ): posts cannot be edited once approved;
content is released under CC-BY; R code goes in plain `<pre>` tags; images
are uploaded as files, not pasted from knitr output.

1. Register at <https://r-posts.com/register/> and log in.
2. Start a new post. Title:
   *DSGE models in R: estimating a New Keynesian model on US data with dsge*
3. Switch the editor to the **Text** tab and paste the whole of
   `wordpress.html`.
4. Switch to **Visual**. At each `[UPLOAD IMAGE: ...]` marker, drag in the
   named PNG from `figures/`, then delete the marker.
5. Preview the post, check the code blocks and images, then submit it for
   review.
