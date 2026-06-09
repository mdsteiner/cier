# Independent reference for the cier_screen() vote combiner.
#
# cier_screen() is an orchestrator, not a new statistic: it runs the index
# battery, lays out a per-index flag table, and then COLLAPSES correlated
# indices to one vote before counting agreement. The only genuinely new logic is
# that collapse + count, so this reference re-derives it from scratch and never
# calls the production combiner (`collapse_votes()` / `cier_screen()`).
#
# The collapse rule: indices sharing a `vote_group` fuse into one vote that fires
# when ANY member flagged the respondent (logical OR); an abstaining member
# (`NA`) contributes FALSE (it did not flag), never TRUE. The per-respondent
# flag count is then the number of vote groups that fired.
#
# The Poisson-binomial agreement itself is already pinned by its own enumeration
# oracle (ref-poisson-binomial-enumeration.R); the screen tests reuse the proven
# `flag_agreement()` on these independently-collapsed votes, so this file owns
# only the collapse + count.

# Collapse a per-index flag table to per-vote-group votes by OR, with NA -> FALSE.
# Args:
#   flags      - matrix/data.frame, respondents x indices, logical (NA allowed),
#                column names are method ids.
#   vote_group - named character map (names = method ids) giving each index's
#                vote group; members sharing a group fuse to one vote.
# Returns: a data.frame, respondents x vote groups (group first-appearance order
#   across the flag columns), logical with no NA.
ref_collapse_votes <- function(flags, vote_group) {
  flags <- as.data.frame(flags, stringsAsFactors = FALSE)
  methods <- colnames(flags)
  groups <- unique(unname(vote_group[methods]))           # first-appearance order
  out <- lapply(groups, function(g) {
    members <- methods[vote_group[methods] == g]
    acc <- rep(FALSE, nrow(flags))
    for (mth in members) {
      f <- flags[[mth]]
      f[is.na(f)] <- FALSE                                 # abstain -> not flagged
      acc <- acc | f                                       # logical OR across members
    }
    acc
  })
  names(out) <- groups
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

# Per-respondent count of vote groups that fired (the agreement count input).
ref_screen_n_flags <- function(votes) {
  as.integer(rowSums(as.data.frame(votes, stringsAsFactors = FALSE)))
}
