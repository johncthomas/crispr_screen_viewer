from crispr_screen_viewer.util import *

import dash_core_components as dcc
import dash_html_components as html
from volcano import spawn_volcanoes
from scatter import spawn_scatter

# this should get a list of valid analyses from a given dir, but for now it'll be hard coded

screens_dir = '/Users/johnc.thomas/Dropbox/crispr/screens_analysis/'

analyses = ['david_756+ddrV2/dav_756-73/take2_libV2_1',
            'harv_dom_RPA_flowscreen/hardom_789-94/take1']

def launcher():
    """User specifies the data set and desired chart"""