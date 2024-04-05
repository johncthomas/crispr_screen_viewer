import unittest
from unittest import TestCase
from glob import glob

from sqlalchemy.orm import Session

from crispr_screen_viewer.database import GeneTable

from crispr_screen_viewer.update_database import *
from crispr_screen_viewer.dataset import ANALYSESTYPES, MetadataTables
from crispr_screen_viewer.database import *
from crispr_screen_viewer.functions_etc import get_resource_path, get_ith_from_all
from crispr_screen_viewer.update_database import create_engine_with_schema, upsert_records
from crispr_screen_viewer.tests.utilities import (
    create_test_database,
    TEST_DB_DIR,
    INFOS,
    INFOS_branched,
    DATA_DIR
)

create_test_database()

class TestDatabaseValues(TestCase):
    """Create database and check some values"""
    def setUp(self):
        engine, metadata = load_test_db_data()

        self.engine = engine
        self.metadata = metadata

    def test_check_stat_tables(self):
        """Checks all FDR values are there and match between source tables and
        database."""
        sample = 'CTRL-TREAT'
        for ans in ANALYSESTYPES:
            df = pd.read_csv(
                get_resource_path(
                    f'tests/test_data/flat_style/results/test1/tables/result.{ans.name}_table.csv',
                ),
                index_col=0, header=[0,1]
            ).sort_index()
            fdrs = df.loc[:, (sample, 'fdr')]

            fdrs = fdrs.loc[~fdrs.index.str.startswith('Non-')]

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
        expected = {
            'ATR-KO',
            'ATR-KO (with Drug)',
            'Drug',
            'Drug (in ATR-KO cells)',
            'No treatment',
            'Named contrast'
        }

        # with Session(self.engine) as S:
        #     observed = S.query(ComparisonTable.treatment_label).distinct().all()
        #     observed = set([t[0] for t in observed])

        observed = set(self.metadata.comparisons['Treatment'].unique())

        self.assertTrue(observed == expected)

    def test_gene_labels(self):

        for libd, gene in [('lib1_brunello', 'ATR',), ('lib2_mouse','Adgrl2')]:

            fn_exo = get_resource_path(f'tests/test_data/exorcise_libs_A/{libd}/exorcise.tsv.gz')
            exo_table = pd.read_csv(fn_exo, sep='\t', index_col=0)

            # Don't want to hardcode the correct values since the tables could change,
            #   so just do a quick check for a matching value
            external_id = exo_table.loc[exo_table.exo_symbol == gene, 'inherit.externalId'].values[0]
            with Session(self.engine) as session:
                res = session.query(
                    GeneTable.symbol_with_ids
                ).where(
                    GeneTable.symbol == gene
                ).one()
            gene_label = res[0]
            self.assertTrue(
                external_id in gene_label,
                f"External ID not found in gene label: {gene=}, {external_id=}, {gene_label=}"
            )


class TestDBManipulations(TestCase):
    def setUp(self):
        engine, metadata = load_test_db_data()
        self.engine:Engine = engine
        self.metadata: MetadataTables = metadata
    def test_remove_exp(self):
        def get_expid(session, metadata):
            xpid_db = set(get_ith_from_all(
                session.query(StatTable.experiment_id).distinct().all()
            ))
            xpid_exp = set(metadata.experiments['Experiment ID'].unique())
            xpid_cmp = set(metadata.comparisons['Experiment ID'].unique())

            return xpid_db, xpid_exp, xpid_cmp


        with Session(self.engine) as session:
            xpids = get_expid(session, self.metadata)

            self.assertTrue(xpids[0] == xpids[1] == xpids[2],
                            f"Experiment IDs in test data don't match before manipulation, {xpids}")

            remove_exp = 'test2'

            remaining_exp = xpids[0]
            remaining_exp.remove(remove_exp)

            metadata2 = remove_experiments(
                exp_ids=[remove_exp],
                metadata_tables=self.metadata,
                session=session,
            )

            xpids2 = get_expid(session, metadata2)
            self.assertTrue(xpids2[0] == xpids2[1] == xpids2[2],
                            f"Experiment IDs in test data don't match after manipulation, {xpids}")

            self.assertFalse(
                any([remove_exp in x for x in xpids2]),
                 f"Removed experiment ID was found in at least one place {xpids2}"
            )

            self.assertTrue(
                all([remaining_exp == xpids2[0]]),
                f"Missing experiment ID after removing {remove_exp} - expected: {remaining_exp}, found {xpids2[0]}"
            )


class TestUpdateExp(TestCase):
    def test(self):
        replace_list = INFOS_branched[-1:] # get last exp, as a list
        replace_expid = replace_list[0].experiment_id
        engine, metadata = load_test_db_data()

        self.assertTrue((metadata.comparisons['Experiment ID'] == replace_expid).any())
        self.assertTrue((metadata.experiments['Experiment ID'] == replace_expid).any())

        update_database(
            TEST_DB_DIR,
            replace_list,
            update_experiments=True,
        )

        engine, metadata = load_test_db_data()

        for md in ('experiments', 'comparisons'):
            xpids = set(getattr(metadata, md)['Experiment ID'].unique())
            self.assertTrue(
                replace_expid in xpids,
                f"Missing experiment ID {replace_expid} in {md} expids: {xpids}, replacing {replace_list}"
            )


class TestUpsert(TestCase):
    def test_upsert(self):
        engine = create_engine_with_schema()
        with Session(engine) as session:
            upsert_records(
                [{'id': 'g1', 'symbol': 'symb1', 'organism':'Human'}],
                session,
                GeneTable
            )

            upsert_records(
                [
                    {'id': 'g1', 'symbol': 'symb1b', 'organism':'Human'},
                    {'id': 'g2', 'symbol': 'symb2', 'organism':'Human' }
                ],
                session,
                GeneTable
            )

            symbols = set(get_ith_from_all(session.query(GeneTable.symbol).all()))

            self.assertTrue(symbols == {'symb1b', 'symb2'})

def test_add_genes_from_symbols():
    from crispr_screen_viewer.update_database import (
        create_engine_with_schema,
        load_stats_csv,
        add_genes_from_symbols,
        get_gene_symbols_in_db
    )
    from importlib import resources
    ngn = create_engine_with_schema(echo=True)

    fn = resources.files("crispr_screen_viewer").joinpath("tests/test_data/branched_style/test1/re/res/tables/result.drugz_table.csv").__str__()
    test_table = load_stats_csv(fn)

    with Session(ngn) as S:
        add_genes_from_symbols(test_table.index[:-5], 'Human', S)
        assert len(get_gene_symbols_in_db(S)) == len(test_table.index) - 5
        add_genes_from_symbols(test_table.index, 'Human', S)
        assert len(get_gene_symbols_in_db(S)) == len(test_table.index)




if __name__ == '__main__':
    unittest.main()



