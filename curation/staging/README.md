# Active Staging

Place YAML records here when you are actively testing or applying them with
`scripts/stage_records.py`.

The staging script only reads files explicitly passed with `--staging`.

Generated staging files should expose the complete template field set for each
record kind, with unknown optional values left blank for curator review.
