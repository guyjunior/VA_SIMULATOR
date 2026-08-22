"""SQLite wearing a MariaDB face.

WHY THIS LAYER EXISTS
---------------------
The analytical pipeline (`app/services/virtual_analyst.py`) is a FAITHFUL PORT
of the code running in production against MariaDB. The less it changes, the
smaller the chance its results diverge from the real system - which is the only
thing that matters here.

So instead of rewriting dozens of queries, this layer makes SQLite speak the
dialect that code already speaks:

    cursor = conn.cursor(dictionary=True)   # rows as dicts
    cursor.execute("... WHERE id = %s", (x,))
    cursor.lastrowid

KNOWN TRAP
----------
Placeholder translation is a textual replace of `%s` with `?`. That would break
a literal containing `%s` - `LIKE '%s%'`, for instance. There is none today (the
pipeline's LIKE clauses use `'%SUSPECT%'`, with an uppercase S, which the
replace leaves alone), and the way out is: pass the literal as a PARAMETER,
never inline in the SQL.
"""

from __future__ import annotations

import datetime as _dt
import re
import sqlite3
from functools import lru_cache
from pathlib import Path

from app import config

# Python 3.12 deprecated the implicit date/datetime adapters. Registering them
# explicitly silences the DeprecationWarning and, more importantly, pins the
# stored format (ISO with a space, same as MariaDB's DATETIME).
sqlite3.register_adapter(_dt.datetime, lambda v: v.isoformat(sep=" ", timespec="seconds"))
sqlite3.register_adapter(_dt.date, lambda v: v.isoformat())


_NOW_RE = re.compile(r"\b(?:NOW|UTC_TIMESTAMP|CURRENT_TIMESTAMP)\s*\(\s*\)", re.IGNORECASE)
_INSERT_IGNORE_RE = re.compile(r"\bINSERT\s+IGNORE\s+INTO\b", re.IGNORECASE)


@lru_cache(maxsize=512)
def translate(sql: str) -> str:
    """Convert the MySQL dialect used by the port into SQLite's."""
    sql = _NOW_RE.sub("CURRENT_TIMESTAMP", sql)
    sql = _INSERT_IGNORE_RE.sub("INSERT OR IGNORE INTO", sql)
    return sql.replace("%s", "?")


class CompatCursor:
    """A cursor with the signature the ported code expects."""

    def __init__(self, cursor: sqlite3.Cursor, dictionary: bool):
        self._cursor = cursor
        self._dictionary = dictionary

    # -- execution -----------------------------------------------------------
    def execute(self, sql, params=None):
        self._cursor.execute(translate(sql), tuple(params) if params else ())
        return self

    def executemany(self, sql, seq):
        self._cursor.executemany(translate(sql), [tuple(p) for p in seq])
        return self

    # -- reading -------------------------------------------------------------
    def _wrap(self, row):
        if row is None:
            return None
        return dict(row) if self._dictionary else tuple(row)

    def fetchone(self):
        return self._wrap(self._cursor.fetchone())

    def fetchall(self):
        return [self._wrap(r) for r in self._cursor.fetchall()]

    def fetchmany(self, size=None):
        rows = self._cursor.fetchmany(size) if size else self._cursor.fetchmany()
        return [self._wrap(r) for r in rows]

    def __iter__(self):
        for row in self._cursor:
            yield self._wrap(row)

    # -- metadata ------------------------------------------------------------
    @property
    def lastrowid(self):
        return self._cursor.lastrowid

    @property
    def rowcount(self):
        return self._cursor.rowcount

    @property
    def description(self):
        return self._cursor.description

    def close(self):
        self._cursor.close()


class CompatConnection:
    """A connection with `cursor(dictionary=...)`, like mysql-connector's."""

    def __init__(self, raw: sqlite3.Connection):
        self._raw = raw

    def cursor(self, dictionary: bool = False, **_ignored) -> CompatCursor:
        return CompatCursor(self._raw.cursor(), dictionary)

    def commit(self):
        self._raw.commit()

    def rollback(self):
        self._raw.rollback()

    def close(self):
        self._raw.close()

    def is_connected(self) -> bool:
        try:
            self._raw.execute("SELECT 1")
            return True
        except sqlite3.Error:
            return False

    @property
    def raw(self) -> sqlite3.Connection:
        return self._raw

    def __enter__(self):
        return self

    def __exit__(self, exc_type, *_):
        if exc_type is None:
            self.commit()
        else:
            self.rollback()
        self.close()
        return False


def connect(path: Path) -> CompatConnection:
    """Open a connection.

    One connection per thread - a batch is processed on a background thread
    while Flask keeps serving, and a SQLite connection does not safely cross
    threads.

    WAL plus busy_timeout is what lets the browser READ the task (status,
    partial results) while the pipeline WRITES to it. Without WAL, every read
    during processing would come back with 'database is locked'.
    """
    conn = sqlite3.connect(str(path), timeout=30.0, isolation_level=None)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = NORMAL")
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA busy_timeout = 30000")
    return CompatConnection(conn)


def lake_conn() -> CompatConnection:
    """LAKE database: ingested samples, channels, scan filters seen in the mzML."""
    return connect(config.LAKE_DB)


def analysis_conn() -> CompatConnection:
    """ANALYSIS database: method, substances, tasks, results."""
    return connect(config.ANALYSIS_DB)
