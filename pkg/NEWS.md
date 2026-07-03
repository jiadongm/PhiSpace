# PhiSpace 1.1.0.9000

* `scoreCells()` now calculates each class's correlation score using that
  class's effective signature genes. Previous versions used every gene shared
  by the reference and query.
* `scoreCells()` gains a `signature_genes` argument for supplying a character
  vector to a one-class reference or an exactly named list to a multi-class
  reference.
* Runs now record whether signatures were generated or supplied and which
  effective genes were used for correlation scoring.
