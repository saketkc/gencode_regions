"""Top-level package for gtfbase."""

__author__ = """Saket Choudhary"""
__email__ = 'saketkc@gmail.com'
__version__ = '0.1.0'

from .ensembl_data_manager import EnsemblDataManager
from .gene_db import GeneDB
from .bed_tool_builder import BedToolBuilderFactory
from .consts import TEMP_DIR_NAME
from .gtfbase_cli import main
