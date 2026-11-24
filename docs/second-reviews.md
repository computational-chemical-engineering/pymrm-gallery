# Second reviews

Every page gets one independent review before it is merged (`gallery-reviewer`, Opus). A
second review by a different or newer model is done **after** the merge, on selected pages,
and its verified must-fix items are fixed in a follow-up PR. It is not a gate for publishing.

Decided by the maintainer on 2026-11-24, after the reviewer test C5: on three new pages
(F1.8, F2.1, J2.6) the Opus and Fable reviews of the same pre-fix commit each found verified
must-fix items the other missed (Fable 17, Opus 13, 24 distinct defects), so a second review adds
quality, but it does not need to hold up publication.

Which pages get one:

1. New pages, starting with J1.7, J2.7 and J5.6.
2. Critical pages published earlier (chosen by the maintainer or the coordinating session).
3. Once the gallery is nearly complete: a slow rolling scan with newer Claude releases.

`second-reviews.csv` records every second review, so it is known which pages have one:

| column | meaning |
|---|---|
| page | catalogue ID |
| model | model ID of the reviewer, as its transcript shows it |
| kind | `post-merge`, or `pre-merge (C5 test)` for the three test pages |
| reviewed_commit | the commit the reviewer read |
| date | date of the review |
| must_fix_verified | must-fix items confirmed by the coordinating session |
| fix_pr | PR that applied the fixes (empty while open; `none` if nothing needed fixing) |
| notes | short context |

A page counts as checked by a model when it has a row with that model and its fix PR is
merged (or `none`).
