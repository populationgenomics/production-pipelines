from cpg_utils import to_path
from cpg_utils.config import prepend_config_paths
from .workflow import (
    get_workflow,
    get_batch,
    get_cohort,
)


defaults_config_path = to_path(__file__).parent / 'defaults.toml'
assert defaults_config_path.exists(), defaults_config_path
prepend_config_paths([str(defaults_config_path)])
