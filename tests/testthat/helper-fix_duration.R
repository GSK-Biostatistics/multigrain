# helper function to be used with the `expect_snapshot()` `transform` argument.
# it takes a string and replaces measured duration (e.g. `34ms`, `1.2s`,
# `1m 18.9s`) with a fixed `"10ms"`. duration is expected to fluctuate and
# we're not interested in the exact number, rather the milestone messaging
fix_duration <- function(out) {
    gsub(
        # `\[` and `\]` match literal square brackets around the duration.
        # `(?:A|B)` means either:
        #   A: `[0-9]+(?:\.[0-9]+)?ms` for milliseconds, e.g. `34ms`, `12.5ms`
        #   B: `(?:[0-9]+m\s*)?[0-9]+(?:\.[0-9]+)?s` for seconds with optional
        #      minutes prefix, e.g. `1.8s`, `18s`, `1m 18.9s`, `1m18.9s`
        pattern = "\\[(?:[0-9]+(?:\\.[0-9]+)?ms|(?:[0-9]+m\\s*)?[0-9]+(?:\\.[0-9]+)?s)\\]", # nolint: line_length_linter
        replacement = "[10ms]",
        x = out,
        perl = TRUE
    )
}
