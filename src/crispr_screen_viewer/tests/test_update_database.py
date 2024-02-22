import unittest
from unittest import TestCase

from crispr_screen_viewer import update_database

class TestAddData(TestCase):
    def setUp(self):
        self.engine = update_database.create_engine_with_schema()

    def test_add_data_exorcise_v1(self):
        infos = update_database.get_paths_simons_structure_v1(['./test_data/exorcise_style/test1'])
        update_database.add_data_to_database(infos, self.engine, overwrite=False)


if __name__ == '__main__':
    unittest.main()