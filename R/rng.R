# Local-seed RNG helpers (with_local_seed): a non-NULL seed leaves the global
# `.Random.seed` byte-identical to its pre-call state; a NULL seed is a transparent
# pass-through drawing from the ambient stream.

# Restore the captured pre-set.seed() state: re-assign saved `.Random.seed`, or remove
# it when none existed.
restore_random_seed <- function(saved) {
  global <- globalenv()
  if (is.null(saved)) {
    if (exists(".Random.seed", envir = global, inherits = FALSE)) {
      rm(".Random.seed", envir = global)
    }
  } else {
    global[[".Random.seed"]] <- saved
  }
}

# Run thunk `fn` under a temporary local seed and return its value. Non-NULL `seed`
# saves, set.seed()s, and restores on exit; NULL `seed` is a pass-through.
with_local_seed <- function(seed, fn) {
  if (!is.null(seed)) {
    saved <- globalenv()[[".Random.seed"]]   # NULL when no RNG has been drawn yet
    on.exit(restore_random_seed(saved), add = TRUE)
    set.seed(seed)
  }
  fn()
}
