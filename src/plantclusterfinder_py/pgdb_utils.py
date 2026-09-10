"""Utils and processing functions for PGDB data."""

import logging
import time
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
        y: dict[str, Any] = {}  # dict.fromkeys(attrib)
        # create a copy of y to start off each row
        # y_prime = y.copy()
        # previous column key, is set when previous key is valid
        prev_key: str | None = None
        for line in fh:
            # data row ends, split lines and reset row
            if line.startswith("//"):
                #  return all attrib values, even if they are not found
                yield {k: y.get(k) for k in attrib}
                y = {}
                continue
            # ignore comment rows
            elif line.startswith("#"):
                continue
            # value continues from previous line
            elif line.startswith("/"):
                if prev_key is not None:
                    if not isinstance(y.get(prev_key), list):
                        y[prev_key] = [y[prev_key], line[1:]]
                        continue
                    last_val = y[prev_key][-1]
                    last_val += line[1:]
                    y[prev_key][-1] = last_val
                else:
                    # prev_key is not set, ignore line
                    continue
            # get key, value pair from line
            elif " - " in line:
                key, *value_arr = line.split(" - ")
                value: str = (" - ".join(value_arr)).strip()
                # value is an annotation-value pair
                if line.startswith("^") and prev_key is not None and prev_key in attrib:
                    attrib_key = prev_key + "_attrib_" + key
                    if attrib_key not in y or y[attrib_key] is None:
                        y[attrib_key] = []
                    y[attrib_key] = y[attrib_key].append({key: value})
                else:
                    prev_key = key
                    # each row should have exactly one unique_key
                    if key in unique_keys:
                        y[key] = value
                    # only store values in attrib
                    elif key in attrib:
                        y[key] = y.get(key, []) + [value]
                    else:
                        prev_key = None
            else:
                raise ValueError("Invalid PGDB File")

    # The PGDB source data doesn't enforce encoding, but Arabidopsis is latin-1
    with pgdb_filename.open(encoding=pgdb_encoding) as f:
        try:
            ret_df = pl.from_dicts(_row_gen(f, attrib, unique_keys), infer_schema_length=None)
            _log.debug("Extracted from %s, %s", pgdb_filename, ret_df.head())
        except pl.exceptions.ComputeError as e:
            schema = set(unique_keys).union(attrib)
            ret_df = pl.from_dicts(_row_gen(f, attrib, unique_keys), schema=list(schema))
            if any(k not in ret_df.columns for k in unique_keys):
                _log.debug("columns: %s, expected: %s", ret_df.columns, unique_keys)
                time.sleep(1)  # wait to flush log
                raise ValueError("Could not parse PGDB source: can't find unique key") from e
    return ret_df


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
