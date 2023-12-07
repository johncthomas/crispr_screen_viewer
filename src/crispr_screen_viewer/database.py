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



__all__ = "GeneTable", "ExperimentTable", 'ComparisonTable', 'StatTable',  "Base"

# function for creating columns
mcol = mapped_column

MInt = orm.Mapped[int]
MFloat = orm.Mapped[float]
MStr = orm.Mapped[str]
# nullable versions
MIntN = orm.Mapped[Optional[int]]
MFloatN = orm.Mapped[Optional[float]]
MStrN = orm.Mapped[Optional[str]]


class Base(orm.DeclarativeBase):
    pass


# going with a simple
# class AnalysisTypeTable(Base):
#     __tablename__ = 'analysis_type'
#
#     id: MInt = mcol(primary_key=True)
#     name: MStr
#     shortname: MStr
#     label: MStrN


class GeneTable(Base):
    __tablename__ = "gene"

    #id: MInt = mcol(primary_key=True)
    name: MStr = mcol(primary_key=True)
    ensembl: MStrN
    ncbi: MStrN
    hgnc: MStrN
    mgi: MStrN


class ExperimentTable(Base):
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


class ComparisonTable(Base):
    __tablename__ = 'comparison'

    #id: MInt = mcol(primary_key=True)
    stringid: MStr = mcol(primary_key=True)
    experiment: MStr = ForeignKey(ExperimentTable.stringid)
    treatment_label: MStr
    timepoint: MStr
    cell: MStr
    ctrl: MStr
    treat: MStr
    ko: MStr
    dose: MStrN
    gi: MStrN
    days_grown: MIntN
    library: MStr
    analyses_bitmask: MInt


class StatTable(Base):
    __tablename__ = 'stat'
    __table_args__ = (
        # sqla.PrimaryKeyConstraint(
        sqla.UniqueConstraint(
            'comparison_id', 'gene_id', 'analysis_type_id'
        ),
    )

    # composite foreign primary keys, defined in table_args
    #  so if they change name, __table_args__ needs to change
    comparison_id: MInt = orm.mapped_column(ForeignKey(ComparisonTable.stringid), primary_key=True)
    #gene_id: MInt = orm.mapped_column(ForeignKey(GeneTable.id), primary_key=True)
    gene_id: MStr = mapped_column(primary_key=True)
    analysis_type_id: MInt = orm.mapped_column(primary_key=True)

    score: MFloatN
    fdr: MFloatN
    fdr10: MFloatN
    pos_p: MFloatN
    neg_p: MFloatN






