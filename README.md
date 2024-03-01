# Documentation

## Installation
`pip install .`

## Dataset
Input data is currently several CSV files and a SQL database. In the future this is likely to move to a single 
database for everything. Example dataset is in the data directory.

### Statistics tables
Wide tables of values for each gene (row) and comparison (column). File names are suffixed with mag|drz
indicating MAGeCK or DrugZ analyses, and suffixed with pos_p|neg_p|fdr|score indicating p-value
of increasing/decreasing abundance, FDR corrected p-value, or "score" – LFC|NormZ for MAGeCK|DrugZ analyses.

These tables are generated by https://github.com/johncthomas/crispr_tools script `pipeline_to_viewer.py` if
using that package.

### Metadata tables
`comparisons_metadata.csv`: Information about comparisons in statistics tables, e.g. treatments and cell types.
`experiments_metadata.csv`: Information about experiments referenced in comparisons metadata, e.g. citation.
`previous_and_id.csv`: Maps verious database IDs and previous symbols used to the symbols used in the statistics
tables.

See the files in the data directory for structure.

### Database
`database.db` is created from the statistics tables using the script `create_db.py`.

## Running
### via script
After installing the package run the following for 
`launch-crsv --help`

### Gunicorn
Create a python module to create the server object, called (for example) `test_crsv.py`
```python
from crispr_screen_viewer import launch  
from importlib import resources 

# locate the included test data
data_dir = resources.files("crispr_screen_viewer").joinpath("data").__str__()  
db_url = f"sqlite:///{data_dir}/database.db"  

app = launch.init_app(  
    data_path=data_dir,  
    database_url=db_url,  
    debug_messages=True  
)  

server = app.server
```

Then  `gunicorn --bind 0.0.0.0:8050 test_crsv:server`

### Gunicorn via docker
Put your data in the `src/crispr_screen_viewer/data` directory before building the image, then run 
```sh
docker run \
	--user appuser \
	-p 8050:8050 \
	--workdir /app/src/crispr_screen_viewer/ \
	<docker image> \
	gunicorn --bind 0.0.0.0:8050 launch:"get_server()"
```