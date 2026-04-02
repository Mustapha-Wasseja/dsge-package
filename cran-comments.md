## Resubmission

Addressing reviewer feedback from Benjamin Altmann:

- Added \value tags to all exported .Rd files
- Replaced \dontrun{} with \donttest{} for long-running examples
- Unwrapped fast examples (< 5 sec) so they run without wrappers
- Fixed incomplete examples to be fully self-contained

## R CMD check results

0 errors | 0 warnings | 0 notes
