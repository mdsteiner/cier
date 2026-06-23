# Tests for the local-seed RNG helper (R/rng.R).
#
# with_local_seed(seed, fn) applies `seed` for the duration of the thunk `fn`
# WITHOUT disturbing the caller's RNG stream (the save / set.seed / on.exit-restore
# idiom shared by kernel_rpr, resolve_personfit_null_cutoff, and cier_simulate); a
# NULL seed is a transparent pass-through that draws from the ambient stream and
# advances it (no save, no reseed). restore_random_seed() is exercised through
# both the saved-then-restored path and the no-prior-seed (remove-on-exit) path.

test_that("a seeded call leaves the ambient .Random.seed bytewise untouched", {
  withr::with_seed(123L, {
    before <- get(".Random.seed", envir = globalenv())
    out <- with_local_seed(42L, function() runif(3L))
    after <- get(".Random.seed", envir = globalenv())
    expect_identical(after, before)            # caller's stream undisturbed
    expect_length(out, 3L)                      # the thunk still ran under the seed
  })
})

test_that("with no prior .Random.seed, a seeded call removes it again on exit", {
  withr::with_preserve_seed({
    if (exists(".Random.seed", envir = globalenv())) {
      rm(".Random.seed", envir = globalenv())
    }
    with_local_seed(7L, function() runif(2L))
    # No seed existed before the call, so restore must return to that absent state
    # (a naive save/restore that always re-assigns would leak a seed here).
    expect_false(exists(".Random.seed", envir = globalenv()))
  })
})

test_that("a NULL seed draws from and advances the ambient stream", {
  withr::with_seed(99L, {
    before <- get(".Random.seed", envir = globalenv())
    val <- with_local_seed(NULL, function() runif(1L))
    after <- get(".Random.seed", envir = globalenv())
    expect_false(identical(before, after))      # no save/restore: stream advanced
    expect_length(val, 1L)
  })
})

test_that("a NULL seed reseeds nothing (reproducible only under an outer seed)", {
  # A hardcoded internal seed would make a == c; a set.seed(NULL) entropy reseed
  # would make a != b. Neither: a NULL seed defers entirely to the ambient stream.
  a <- withr::with_seed(5L, with_local_seed(NULL, function() runif(4L)))
  b <- withr::with_seed(5L, with_local_seed(NULL, function() runif(4L)))
  cc <- withr::with_seed(6L, with_local_seed(NULL, function() runif(4L)))
  expect_identical(a, b)
  expect_false(isTRUE(all.equal(a, cc)))
})

test_that("a seeded call is reproducible regardless of the ambient stream", {
  withr::with_preserve_seed({
    a <- with_local_seed(2026L, function() runif(5L))
    invisible(runif(10L))                        # perturb the ambient stream
    b <- with_local_seed(2026L, function() runif(5L))
    expect_identical(a, b)                        # seed makes the draw deterministic
  })
})

test_that("a seeded call applies exactly set.seed(seed), nothing more", {
  # Pins the draw to set.seed(seed) verbatim, killing any consistent seed
  # mistransform (set.seed(seed + 1L), as.integer(), abs(), ...) that a
  # reproducibility-only test cannot see.
  ref <- withr::with_preserve_seed({
    set.seed(7L)
    runif(3L)
  })
  expect_identical(with_local_seed(7L, function() runif(3L)), ref)
})

test_that("with_local_seed returns the thunk's value (both seed modes)", {
  expect_identical(with_local_seed(1L, function() "ok"), "ok")
  expect_identical(with_local_seed(NULL, function() 42L), 42L)
})
