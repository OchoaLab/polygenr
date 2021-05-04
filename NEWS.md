# 2021-05-04 - polygenr 0.0.0.9000

- Initial commit, includes functions `glmnet_pca`, `scores_glmnet`, `anova2`, and `anova_glmnet`.
  Functions are well documented.

# 2021-05-04 - polygenr 0.0.1.9000

- Functions `anova2` and `anova_glmnet` now work when genotype matrix `X` has rownames with special characters.
  Before, for example if a locus name had the common format "chr:pos" where chr and pos are both integers, these functions would fail with this error: "Error in terms.formula(formula, data = data) : invalid model formula in ExtractVars".
  Now all variables are quoted with backticks, which avoids this problem.
  Tested to work when names have math operations such as "a+b" (treated as a single variable name now).
  Additionally, an informative error is thrown if a locus name actually contains a backtick (only special character we can't handle).
