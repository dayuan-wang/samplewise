# docs/ci/

`R-CMD-check.yaml` is the GitHub Actions workflow that runs `R CMD check` and
the tests on every push. It lives here instead of `.github/workflows/` because
the token used for the first push (2026-09-27) lacked the `workflow` scope.
To activate it:

```bash
gh auth refresh -h github.com -s workflow
git mv docs/ci/R-CMD-check.yaml .github/workflows/R-CMD-check.yaml
git commit -m "Enable CI" && git push
```
