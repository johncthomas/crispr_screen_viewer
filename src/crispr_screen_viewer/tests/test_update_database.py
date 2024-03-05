import os
import unittest, math
from unittest import TestCase

import pandas as pd
import numpy as np

from sqlalchemy.orm import Session

from crispr_screen_viewer import update_database
from crispr_screen_viewer.dataset import ANALYSESTYPES, DataSet
from crispr_screen_viewer.database import *


import logging
logger = logging.getLogger(__name__)
logging.basicConfig()

logging.getLogger('update_database').setLevel('INFO')


class TestDatabaseExorciseV1(TestCase):
    """Tests using input in Exorcise output structure."""
    def setUp(self):

        outdir = './test_data/exorcise_out'

        import glob

        files = glob.glob(outdir+'/*')
        for f in files:
            os.remove(f)

        self.outdir = outdir
        self.engine = update_database.create_engine_with_schema()
        infos = update_database.get_paths_simons_structure_v1(['./test_data/exorcise_style/test1'])
        update_database.add_data_to_database(
            infos,
            self.engine,
            outdir=outdir,
            overwrite=False
        )
        # self.comparisons_metadata = pd.read_csv(os.path.join(outdir, 'comparisons_metadata.csv.gz'), index_col=)
        # self.experiments_metadata = pd.read_csv(os.path.join(outdir, 'experiments_metadata.csv.gz'))

    def test_check_stat_tables(self):
        """Checks all FDR values are there and match between source tables and
        database."""
        sample = 'CTRL-TREAT'
        for ans in ANALYSESTYPES:
            df = pd.read_csv(
                f'./test_data/exorcise_style/test1/re/res/tables/result.{ans.name}_table.csv',
                index_col=0, header=[0,1]
            ).sort_index()
            fdrs = df.loc[:, (sample, 'fdr')]

            with Session(self.engine) as S:
                db_pvals = S.query(StatTable.fdr).where(
                    StatTable.analysis_type_id == ans.id,
                    StatTable.comparison_id == f'test1.{sample}',
                ).order_by(
                    StatTable.gene_id
                ).all()

                self.assertTrue(all(
                        np.isclose(
                            np.array(db_pvals).ravel(),
                            fdrs.values
                        )
                ))

    def test_treatment_labels(self):
        expected = {'ATR-KO',
         'ATR-KO (with Drug)',
         'Drug',
         'Drug (in ATR-KO cells)',
         'No treatment'}

        with Session(self.engine) as S:
            observed = S.query(ComparisonTable.treatment_label).distinct().all()
            observed = set([t[0] for t in observed])

        self.assertTrue(observed == expected)


def test_add_genes_from_symbols():
    from crispr_screen_viewer.update_database import (
        create_engine_with_schema,
        load_stats_csv,
        add_genes_from_symbols,
        get_gene_symbols_db
    )
    from importlib import resources
    ngn = create_engine_with_schema(echo=True)

    fn = resources.files("crispr_screen_viewer").joinpath("tests/test_data/exorcise_style/test1/re/res/tables/result.drugz_table.csv").__str__()
    test_table = load_stats_csv(fn)

    with Session(ngn) as S:
        add_genes_from_symbols(test_table.index[:-5], 'Human', S)
        assert len(get_gene_symbols_db(S)) == len(test_table.index)-5
        add_genes_from_symbols(test_table.index, 'Human', S)
        assert len(get_gene_symbols_db(S)) == len(test_table.index)




if __name__ == '__main__':
    unittest.main()
