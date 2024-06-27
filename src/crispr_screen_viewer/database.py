import dataclasses
import typing
from typing import (
    Mapping,
    Optional,
)
from datetime import datetime

import sqlalchemy as sqla
from sqlalchemy import (
    ForeignKey, ForeignKeyConstraint, select, insert, Engine
)
from sqlalchemy import orm
from sqlalchemy.orm import Session, mapped_column, Mapped



__all__ = "ExperimentTable", 'ComparisonTable', 'StatTable', "TableBase", 'GeneTable'

# function for creating columns
mcol = mapped_column

# declaritive type shorthands
MInt = orm.Mapped[int]
MFloat = orm.Mapped[float]
MStr = orm.Mapped[str]
# nullable types
MIntN = orm.Mapped[Optional[int]]
MFloatN = orm.Mapped[Optional[float]]
MStrN = orm.Mapped[Optional[str]]


class TableBase(orm.DeclarativeBase):
    pass


class GeneTable(TableBase):
    __tablename__ = "gene"

    id: MStr = mcol(primary_key=True)
    symbol: MStr = mcol(unique=True)
    official_id: MStrN
    symbol_with_ids: MStrN
    organism: MStr


class ExperimentTable(TableBase):
    __tablename__ = 'experiment'

    #id: MInt = mcol(primary_key=True)
    stringid: MStr = mcol(primary_key=True)
    date: orm.Mapped[Optional[datetime]]
    library: MStr
    doi: MStrN
    representation: MStrN
    moi: MStrN
    description: MStrN
    notes: MStrN
    reference: MStrN
    source: MStrN
    citation: MStrN


class ComparisonTable(TableBase):
    __tablename__ = 'comparison'

    #id: MInt = mcol(primary_key=True)
    stringid: MStr = mcol(primary_key=True)
    experiment: MStr = ForeignKey(ExperimentTable.stringid)
    contrast: MStr
    treatment_label: MStr
    timepoint: MStr
    cell: MStr
    control_sample: MStr
    test_sample: MStr
    ko: MStr
    control_treatment: MStrN
    control_ko: MStrN
    dose: MStrN
    gi: MStrN
    days_grown: MIntN
    library: MStr
    notes:MStrN



class StatTable(TableBase):
    __tablename__ = 'stat'
    __table_args__ = (
        sqla.UniqueConstraint(
            'comparison_id', 'gene_id', 'analysis_type_id'
        ),
    )

    # composite foreign primary keys, defined in table_args
    #  so if they change name, __table_args__ needs to change
    comparison_id: MInt = mcol(ForeignKey(ComparisonTable.stringid), primary_key=True)
    experiment_id: MStr = mcol(ForeignKey(ExperimentTable.stringid))
    gene_id: MInt = orm.mapped_column(ForeignKey(GeneTable.id), primary_key=True)
    analysis_type_id: MInt = mcol(primary_key=True)

    score: MFloatN
    fdr: MFloatN
    fdr10: MFloatN
    pos_p: MFloatN
    neg_p: MFloatN






