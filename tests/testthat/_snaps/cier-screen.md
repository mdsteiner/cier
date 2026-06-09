# print renders the locked cli summary (default battery)

    Code
      print(sc)
    Output
      -- cier_screen -----------------------------------------------------------------
      8 indices on 300 respondents; 7 votes after collapsing shared constructs.
      
      Per-index flags (transparent; not independent votes):
        cier_longstring            0 / 300 (0.0%)
        cier_irv                   15 / 300 (5.0%)
        cier_even_odd              45 / 300 (15.0%)  [consistency]
        cier_personal_reliability  15 / 300 (5.0%)  [consistency]
        cier_psychsyn              0 / 0 (--%)
        cier_psychant              0 / 0 (--%)
        cier_mahalanobis           0 / 300 (0.0%)
        cier_person_total          15 / 300 (5.0%)
      
      Cross-index agreement (observed vs independence baseline):
        flagged by >= 1 vote: 74 / 300 (24.7%); expected 24.8%
        flagged by >= 2 votes: 5 / 300 (1.7%); expected 1.8%
        flagged by >= 3 votes: 1 / 300 (0.3%); expected 0.0%  <- excess
        flagged by >= 4 votes: 0 / 300 (0.0%); expected 0.0%
        flagged by >= 5 votes: 0 / 300 (0.0%); expected 0.0%
        flagged by >= 6 votes: 0 / 300 (0.0%); expected 0.0%
        flagged by >= 7 votes: 0 / 300 (0.0%); expected 0.0%
      
      Skipped: 0
      

# print reports skipped methods with their reasons

    Code
      print(sc)
    Output
      -- cier_screen -----------------------------------------------------------------
      6 indices on 60 respondents; 6 votes after collapsing shared constructs.
      
      Per-index flags (transparent; not independent votes):
        cier_longstring            0 / 60 (0.0%)
        cier_irv                   3 / 60 (5.0%)
        cier_psychsyn              0 / 0 (--%)
        cier_psychant              0 / 0 (--%)
        cier_mahalanobis           0 / 60 (0.0%)
        cier_person_total          3 / 60 (5.0%)
      
      Cross-index agreement (observed vs independence baseline):
        flagged by >= 1 vote: 6 / 60 (10.0%); expected 9.8%  <- excess
        flagged by >= 2 votes: 0 / 60 (0.0%); expected 0.2%
        flagged by >= 3 votes: 0 / 60 (0.0%); expected 0.0%
        flagged by >= 4 votes: 0 / 60 (0.0%); expected 0.0%
        flagged by >= 5 votes: 0 / 60 (0.0%); expected 0.0%
        flagged by >= 6 votes: 0 / 60 (0.0%); expected 0.0%
      
      Skipped: 4
        cier_even_odd: needs item metadata (items = NULL)
        cier_personal_reliability: needs item metadata (items = NULL)
        cier_gnormed: needs item metadata (items = NULL)
        cier_ht: needs item metadata (items = NULL)
      

