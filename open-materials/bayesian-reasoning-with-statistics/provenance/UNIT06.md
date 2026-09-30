# Unit 6 fictional workshop data provenance

This 80-row CSV is freshly generated for this book by `scripts/generate_unit06.py` with Python's standard-library `random.Random(20260928)`, eight declared group sizes, and scores rounded to one decimal. Re-running the generator creates identical bytes; the independent verifier does so twice. The generator contains author-side latent centers and a noise scale only to reproduce teaching records. Those inputs are not inferential answers, priors, or known truths about any real students.

Each row is one fictional volunteer's single post-workshop self-report; `workshop_id` records shared setting, not randomized assignment. No human participants, jamovi teaching rows, or external campus data were imported. There is no before measure, missing cell, representative sampling frame, or grounds for workshop effectiveness claims. The resource checksum manifest records this file and its companion scripts.
