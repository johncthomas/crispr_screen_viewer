# Documentation

## Installation
`pip install .`

## Usage
Main script available as `crispr-screen-viewer` with subcommands.

## Manage database
See `crispr-screen-viewer database -h`

Input data is currently several CSV files and a SQL database. In the future this is likely to move to a single 
database for everything. Example dataset is in the data directory.

## Run the viewer
### via script
See `crispr-screen-viewer launch -h`. Requires a database constructed by the above command.

### Gunicorn
Create a python module to create the server object, called (for example) `test_crsv.py`
```python
from crispr_screen_viewer import launch
from crispr_screen_viewer.functions_etc import get_resource_path

# Use the provided example data
db_path = get_resource_path('data/test_db')
app = launch.get_server(data_path=db_path)

server = app
```

Then call `gunicorn --bind 0.0.0.0:8050 test_crsv:server`

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