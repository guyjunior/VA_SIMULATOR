"""ClickHouse clients - one per logical database.

`lake_ch()` reads and writes the raw mzML data; `analysis_ch()` stores the XICs
produced by each run. Same server, different databases.

CONNECTION REUSE: the pipeline calls the insert hundreds of times per sample
(N substances x 4 traces). Opening a client per call exhausts ephemeral ports on
Windows (WinError 10048/10055) - which is why the orchestrator opens ONE client
and passes it down.
"""

from __future__ import annotations

import clickhouse_connect

from app import config


def _client(database: str | None):
    return clickhouse_connect.get_client(
        host=config.CH_HOST,
        port=config.CH_PORT,
        username=config.CH_USER,
        password=config.CH_PASSWORD,
        database=database,
        # The spectra_points insert ships hundreds of thousands of rows per
        # call; the 300s default is short when the machine is busy.
        send_receive_timeout=1800,
    )


def server_ch():
    """Client with NO database attached - required to run CREATE DATABASE."""
    return _client(None)


def lake_ch():
    return _client(config.CH_LAKE_DB)


def analysis_ch():
    return _client(config.CH_ANALYSIS_DB)


def insert_retry(ch, table: str, rows, column_names, max_attempts: int = 5):
    """Insert with backoff.

    It exists for one reason only: MEMORY_LIMIT_EXCEEDED (code 241). A large
    spectra_points batch can hit the server's memory ceiling while it merges;
    waiting and retrying fixes it, because the merge finishes. Any other error
    is NOT retried - it is raised right away.
    """
    import time

    if not rows:
        return

    delay = 5
    for attempt in range(1, max_attempts + 1):
        try:
            ch.insert(table, rows, column_names=column_names)
            return
        except Exception as exc:  # noqa: BLE001 - reclassified below
            text = str(exc)
            transient = "MEMORY_LIMIT_EXCEEDED" in text or "Code: 241" in text
            if not transient or attempt == max_attempts:
                raise
            time.sleep(delay)
            delay = min(delay * 2, 120)
