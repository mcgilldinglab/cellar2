import os

import BTrees.OOBTree
import dash
import dash_bootstrap_components as dbc
import ZODB
import ZODB.FileStorage

from log import setup_logger
from controller.cellar.utils.colors import PALETTE


# In-memory database. Used for storing AnnData objects
db = ZODB.DB(None)
connection = db.open()
dbroot = connection.root
dbroot.adatas = BTrees.OOBTree.BTree()
dbroot.MULTIPLEXER_OUTPUTS = BTrees.OOBTree.BTree()
dbroot.notifications = BTrees.OOBTree.BTree()
dbroot.palettes = BTrees.OOBTree.BTree()
# Using the same palette for every plot right now
# dbroot.palettes['main'] = PALETTE.copy()
# dbroot.palettes['side'] = PALETTE.copy()
dbroot.palettes["universal"] = PALETTE.copy()

# Bootstrap theme
external_stylesheets = [
    dbc.themes.SANDSTONE,
    {
        'href': 'https://use.fontawesome.com/releases/v5.8.1/css/all.css',
        'rel': 'stylesheet',
        'integrity': 'sha384-50oBUHEmvpQ+1lW4y57PTFmhCaXp0ML' +
        '5d60M1M7uH2+nqUivzIebhndOJK28anvf',
        'crossorigin': 'anonymous'
    }
]

app = dash.Dash(
    __name__,
    external_stylesheets=external_stylesheets,
    title="Cellar2",
    suppress_callback_exceptions=True,
    routes_pathname_prefix='/cellar2/',
    requests_pathname_prefix='/cellar2/',
    )
logger = setup_logger('Cellar')
