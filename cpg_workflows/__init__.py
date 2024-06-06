from cpg_utils import to_path
from cpg_utils.config import prepend_config_paths
from cpg_utils.hail_batch import get_batch

from .workflow import (
    get_multicohort,
    get_workflow,
)

defaults_config_path = to_path(__file__).parent / 'defaults.toml'
if defaults_config_path.exists():
    prepend_config_paths([str(defaults_config_path)])
