head -6 LOGS_MS/history.data > history.data
tail -n +7 -q LOGS_*/history.data >> history.data
