"""Utils and processing functions for PGDB data."""

import logging
from collections.abc import Iterator
from pathlib import Path
from typing import IO, Any

import polars as pl

_log = logging.getLogger(__name__)


def read_in_pgdb_flat(
    pgdb_filename: Path,
    attrib: tuple[str, ...],
    unique_keys: tuple[str, ...] = ("UNIQUE-ID",),
    pgdb_encoding: str = "latin-1",
) -> pl.DataFrame:
    """---------------------------------------------------------------------------

      Function that reads in a pgdb flat file (e.g. genes.dat, proteins.dat,
      enzrxns.dat, reactions.dat) saving Attributes given by attrib.

    ---------------------------------------------------------------------------

    Args:
        pgdb_filename: PGDB Flat file that should be read in.
        attrib: Attributes of the PGDB Flat file that should be saved (all
                 others are discarded.)
        unique_keys: Tuple of columns that are unique, i.e. should not be stored as a list.
        pgdb_encoding: Encoding of PGDB file (Default 'latin-1').

    Returns:
        **DataFrame** structure with the contents of each attribute as a column

    """

    # By encapsulating the parsing logic in a generator,
    #  it opens up the door to potentially deferring reading
    #  or simplifying asynchronous reading
    def _row_gen(
        fh: IO, attrib: tuple[str, ...], unique_keys: tuple[str, ...]
    ) -> Iterator[dict[str, Any]]:
        """Generate a dict row from the pgdb flat file."""
        y: dict[str, Any] = dict.fromkeys(attrib)
        prev_key = "COMMENT"
        for line in fh:
            # ignore comment rows, columns
            if line.startswith("#") or line.startswith("COMMENT"):
                if line.startswith("COMMENT - "):
                    prev_key = "COMMENT"
                continue
            # line continues or data row ends
            elif line.startswith("/"):
                # split lines and reset row
                if line.startswith("//"):
                    yield y
                    y = dict.fromkeys(attrib)
                    prev_key = "COMMENT"
                    continue
                # row ends, extend value from previous row
                elif prev_key == "COMMENT":
                    continue
                else:
                    if prev_key not in y:
                        y[prev_key] = None
                        continue
                    elif not isinstance(y[prev_key], list):
                        y[prev_key] += line[1:]
                        continue
                    last_val = y[prev_key][-1]
                    last_val += line[1:]
                    y[prev_key][-1] = last_val

            # get key, value pair from line
            elif " - " in line:
                key, *value_arr = line.split(" - ")
                prev_key = key
                value: str = (" - ".join(value_arr)).strip()
                if key in y:
                    if y[key] is None:
                        if key == unique_keys[0]:
                            y[key] = value
                        else:
                            y[key] = [value]
                    elif not isinstance(y[key], list):
                        y[key] = [y[key], value]
                    else:
                        y[key].append(value)

    # The PGDB source data doesn't enforce encoding, but Arabidopsis is latin-1
    with pgdb_filename.open(encoding=pgdb_encoding) as f:
        return pl.from_dicts(_row_gen(f, attrib, unique_keys), infer_schema_length=None)


def match_pgdb_flat(
    flat1: pl.DataFrame,
    flat2: pl.DataFrame,
    join_on: tuple[str, str],
    primary_id: str = "UNIQUE-ID",
    secondary_id: str = "UNIQUE-ID_right",
) -> pl.DataFrame:
    """---------------------------------------------------------------------------

      Function that maps two flat files of a pgdb to each other.

    ---------------------------------------------------------------------------

    Args:
        flat1: DataFrame of a read out pgdb Flat file.
        flat2: DataFrame of a read out pgdb Flat file.
        join_on: Definition of which attributes in structure 1 should be
                    matched to which attributes of structure 2
        primary_id: Output column for flat1 primary key.
            (Default UNIQUE-ID)
        secondary_id: Output column for flat2 primary key.
            (Default is UNIQUE-ID_right)

    Returns:
        F1_F2: A mapping between items of the first pgdb file to the second
                one.

    """
    _log.info("Matching PGDB flat files %s and %s", flat1, flat2)
    if join_on[0] not in flat1.columns:
        _log.error(
            f"Matching column {join_on[0]} not in input data"
            + f" {flat1.columns}, note the other table's columns are {flat2.columns}"
        )
        raise Exception("Column-matching Mismatch")
    if join_on[1] not in flat2.columns:
        _log.error(
            f"Matching column {join_on[1]} not in input data"
            + f" {flat2.columns}, note the other table's columns are {flat1.columns}"
        )
        raise Exception("Column-matching Mismatch")
    if flat1[join_on[0]].dtype == pl.Null:
        _log.error("Primary table join ID %s is null dtype, cannot join", join_on[0])
        _log.info("Columns of flats are are %s, %s", flat1.columns, flat2.columns)
        raise Exception("Column-matching Mismatch")
    if flat2[join_on[1]].dtype == pl.Null:
        _log.warning(
            "Cannot map %s to %s: dtype is Null, returning Null for second column",
            join_on[0],
            join_on[1],
        )
        _log.info("Columns of flats are are %s, %s", flat1.columns, flat2.columns)
        return flat1.select(pl.col(join_on[0]).alias(primary_id), pl.lit(None).alias(secondary_id))
    flat2 = flat2.explode(join_on[1]) if flat2[join_on[1]].dtype == pl.List else flat2
    flat1 = flat1.explode(join_on[0]) if flat1[join_on[0]].dtype == pl.List else flat1
    # NOTE: Validation on join can prevent data issues, consider adding
    joined_dbs = flat1.join(
        flat2,
        left_on=join_on[0],
        right_on=join_on[1],
        suffix="_right",
    )
    _log.info("Successful join between two PGDB flat file DBs")
    _log.debug("Resulting dataframe has %d rows", joined_dbs.height)
    _log.debug("Mapping UNIQUE-ID (%s) to %s (%s)", primary_id, join_on[0], secondary_id)
    _log.debug("Resulting dataframe: %s", joined_dbs)
    # non_dupe_cols = [x for x in joined_dbs.columns if not x.endswith("_right")]
    return joined_dbs.select(
        pl.col("UNIQUE-ID").alias(primary_id),
        pl.col(join_on[0]).alias(secondary_id),
        #    pl.col(non_dupe_cols),
    )  # .select( pl.col("UNIQUE-ID").alias(primary_id), pl.col(join_on[0]).alias(secondary_id))
