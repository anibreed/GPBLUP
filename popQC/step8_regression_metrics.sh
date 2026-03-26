#!/usr/bin/env bash
set -euo pipefail

# Extract Step-8 regression metrics from popQC log.
# Usage:
#   ./step8_regression_metrics.sh <current_log> [baseline_log]

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <current_log> [baseline_log]" >&2
  exit 1
fi

CUR_LOG="$1"
BASE_LOG="${2:-}"

if [[ ! -f "$CUR_LOG" ]]; then
  echo "ERROR: current log not found: $CUR_LOG" >&2
  exit 1
fi
if [[ -n "$BASE_LOG" && ! -f "$BASE_LOG" ]]; then
  echo "ERROR: baseline log not found: $BASE_LOG" >&2
  exit 1
fi

extract_metrics() {
  local log_file="$1"

  awk '
    function trim(s) { gsub(/^[ \t]+|[ \t]+$/, "", s); return s }
    function lastnum(line,    t) {
      t = line
      gsub(/[^0-9]+/, " ", t)
      t = trim(t)
      split(t, a, /[ ]+/)
      return a[length(a)] + 0
    }

    /Parent update total:/            { pu = lastnum($0) }
    /Multi resolved total:/           { mr = lastnum($0) }
    /ID swap total:/                  { isv = lastnum($0) }
    /Not found total:/                { nf = lastnum($0) }

    /- Changed \/ No change:/ {
      t = $0
      gsub(/[^0-9]+/, " ", t)
      t = trim(t)
      split(t, a, /[ ]+/)
      if (a[1] != "" && a[2] != "") {
        if (chg_lines == 0) { pu_chg = a[1] + 0; pu_nch = a[2] + 0 }
        else if (chg_lines == 1) { mr_chg = a[1] + 0; mr_nch = a[2] + 0 }
        else if (chg_lines == 2) { is_chg = a[1] + 0; is_nch = a[2] + 0 }
        else if (chg_lines == 3) { nf_chg = a[1] + 0; nf_nch = a[2] + 0 }
        chg_lines++
      }
    }

    /MDL=PRE_FAIL->POST_PASS/ { post_pass += lastnum($0) }
    /MDL=PRE_FAIL->POST_FAIL/ { post_fail += lastnum($0) }
    /MDL=PRE_NA/              { pre_na += lastnum($0) }

    END {
      total_chg = pu_chg + mr_chg + is_chg + nf_chg
      total_nch = pu_nch + mr_nch + is_nch + nf_nch
      print "PU=" pu
      print "MR=" mr
      print "IS=" isv
      print "NF=" nf
      print "PU_CHG=" pu_chg
      print "PU_NCH=" pu_nch
      print "MR_CHG=" mr_chg
      print "MR_NCH=" mr_nch
      print "IS_CHG=" is_chg
      print "IS_NCH=" is_nch
      print "NF_CHG=" nf_chg
      print "NF_NCH=" nf_nch
      print "TOTAL_CHG=" total_chg
      print "TOTAL_NCH=" total_nch
      print "POST_PASS=" post_pass
      print "POST_FAIL=" post_fail
      print "PRE_NA=" pre_na
    }
  ' "$log_file"
}

print_block() {
  local title="$1"
  local data="$2"

  echo "[$title]"
  echo "$data" | sed 's/^/  /'
  echo
}

to_assoc() {
  local input="$1"
  while IFS='=' read -r k v; do
    [[ -z "$k" ]] && continue
    eval "$2[\"$k\"]=$v"
  done <<< "$input"
}

CUR_METRICS="$(extract_metrics "$CUR_LOG")"
print_block "Current" "$CUR_METRICS"

if [[ -n "$BASE_LOG" ]]; then
  BASE_METRICS="$(extract_metrics "$BASE_LOG")"
  print_block "Baseline" "$BASE_METRICS"

  declare -A cur
  declare -A base
  to_assoc "$CUR_METRICS" cur
  to_assoc "$BASE_METRICS" base

  echo "[Delta: current - baseline]"
  keys=(PU MR IS NF PU_CHG PU_NCH MR_CHG MR_NCH IS_CHG IS_NCH NF_CHG NF_NCH TOTAL_CHG TOTAL_NCH POST_PASS POST_FAIL PRE_NA)
  for k in "${keys[@]}"; do
    c="${cur[$k]:-0}"
    b="${base[$k]:-0}"
    d=$((c - b))
    printf '  %s=%+d\n' "$k" "$d"
  done
  echo
fi
