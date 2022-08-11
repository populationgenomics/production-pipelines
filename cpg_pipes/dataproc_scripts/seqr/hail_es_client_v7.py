import datetime
import inspect
import logging
import time
from pprint import pformat


try:
    import elasticsearch
except ImportError:
    import os

    os.system("pip install elasticsearch==7.9.1")
    import elasticsearch


handlers = set(logging.root.handlers)
logging.root.handlers = list(handlers)
logger = logging.getLogger()
