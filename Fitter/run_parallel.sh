#!/bin/bash

##############################################
# CONFIGURAZIONE
##############################################

EXEC=./CHeTFitter
CONFIG=macros/edm_tracker.mac

TOTAL_EVENTS=100     # totale degli eventi da processare
NJOBS=2               # numero di job paralleli
OUTDIR_OUT=log/out         # cartella log/output
OUTDIR_ERR=log/err         # cartella log/error

##############################################

mkdir -p "$OUTDIR_OUT"
mkdir -p "$OUTDIR_ERR"

# eventi base per job
EVENTS_PER_JOB=$(( TOTAL_EVENTS / NJOBS ))
REMAINDER=$(( TOTAL_EVENTS % NJOBS ))

echo "===== JOB DISPATCHER ====="
echo "Total events:     $TOTAL_EVENTS"
echo "Jobs:             $NJOBS"
echo "Events/job:       $EVENTS_PER_JOB (+1 for first $REMAINDER jobs)"
echo

for (( job=0; job<NJOBS; job++ )); do

    # alcuni job ottengono 1 evento extra
    extra=0
    if (( job < REMAINDER )); then
        extra=1
    fi

    # calcolo start
    start=$(( job * EVENTS_PER_JOB + (job < REMAINDER ? job : REMAINDER) ))
    end=$(( start + EVENTS_PER_JOB + extra - 1 ))

    # limiti di sicurezza
    if (( start >= TOTAL_EVENTS )); then
        break
    fi
    if (( end >= TOTAL_EVENTS )); then
        end=$(( TOTAL_EVENTS - 1 ))
    fi

    echo "Job $job: processing events $start → $end"

    # LANCIO DEL JOB
    # Sovrascrive i valori nel .mac tramite CLI11
    $EXEC --config "$CONFIG" -q --rangeLoop $start $end \
        > "$OUTDIR_OUT/job_${job}.out" \
        2> "$OUTDIR_ERR/job_${job}.err" &

done

wait
echo
echo "All jobs terminated."
