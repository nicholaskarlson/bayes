# Unit 4 data provenance

`data/unit04_reflection_attention_paired.csv` is a newly generated, entirely fictional 64-row teaching file. `scripts/generate_unit04.py` uses Python's `random.Random(20260927)` and formats one-decimal scores with LF line endings. Its code specifies simulation features for audit; they are **not** known to a reader analyzing the CSV and are not evidence about people. The two designated unusual changes remain valid rows and are not deleted because they complicate a normal likelihood.

One row is one fictional participant, measured at day 0 and day 14. All 64 pairs are complete. There was no random assignment, parallel comparison group, blinded measurement, recruitment frame or construct validation. The design can illustrate a paired mean and conditional posterior; it cannot attribute change to reflection practice or generalize to a population of real students.

The dataset is original to this book. It is not a copy of the jamovi Study 02 records. Its original companion data are CC BY 4.0; `scripts/generate_unit04.py` is included to reproduce the fictional rows. Generator settings disclose construction choices, not truth about real people.
