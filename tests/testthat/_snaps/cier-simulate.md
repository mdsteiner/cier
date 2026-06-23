# print renders the locked summary (checks present)

    Code
      print(s)
    Output
      -- cier_simulate ---------------------------------------------------------------
      Respondents: 200 x 12 items (3 scales) -- 40 careless (20.0%).
      Patterns: alternating 5, diagonal 7, extreme 4, markov 4, midpoint 7, random 4,
      speeder 1, straightline 8.
      Extent: 28 full, 8 partial, 4 temporary.
      Timing: $seconds (totals) + $page_seconds (3 pages: 4+4+4 items).
      Checks: 2 attention checks in $checks (pass sets in $pass).
      Truth: $truth -- careless, pattern, extent, onset_item, offset_item, speeded,
      params.
      Simulated data (power analysis / method comparison / recovery tests), not
      evidence of real-world validity.
      

# print renders the no-checks, all-full variant

    Code
      print(s)
    Output
      -- cier_simulate ---------------------------------------------------------------
      Respondents: 50 x 10 items (2 scales) -- 10 careless (20.0%).
      Patterns: alternating 1, diagonal 1, extreme 2, markov 1, midpoint 2, speeder
      2, straightline 1.
      Extent: 10 full.
      Timing: $seconds (totals) + $page_seconds (2 pages: 5+5 items).
      Checks: none.
      Truth: $truth -- careless, pattern, extent, onset_item, offset_item, speeded,
      params.
      Simulated data (power analysis / method comparison / recovery tests), not
      evidence of real-world validity.
      

# print renders the zero-careless variant

    Code
      print(s)
    Output
      -- cier_simulate ---------------------------------------------------------------
      Respondents: 30 x 10 items (2 scales) -- 0 careless (0.0%).
      Patterns: none.
      Timing: $seconds (totals) + $page_seconds (2 pages: 5+5 items).
      Checks: none.
      Truth: $truth -- careless, pattern, extent, onset_item, offset_item, speeded,
      params.
      Simulated data (power analysis / method comparison / recovery tests), not
      evidence of real-world validity.
      

