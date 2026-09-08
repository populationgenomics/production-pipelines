"""
Tests for restricting the tools that `CramQC` runs via `cramqc/tools`.
"""

from pathlib import Path

import pytest

from cpg_workflows.stages.cram_qc import qc_functions

from .. import set_config
from ..factories.config import PipelineConfig, WorkflowConfig


def config(sequencing_type: str = 'genome', tools: list[str] | None = None) -> PipelineConfig:
    cramqc: dict[str, object] = {'num_pcs': 4}
    if tools is not None:
        cramqc['tools'] = tools

    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='cram-qc-test',
            access_level='test',
            sequencing_type=sequencing_type,  # type: ignore[arg-type]
        ),
        other={'cramqc': cramqc},
    )


def test_all_tools_run_when_unset(tmp_path: Path):
    set_config(config(), tmp_path / 'config.toml')

    assert [qc.name for qc in qc_functions()] == [
        'somalier',
        'verifybamid',
        'samtools_stats',
        'picard_collect_metrics',
        'picard_wgs_metrics',
    ]


def test_restricts_to_the_named_tools(tmp_path: Path):
    set_config(config(tools=['somalier']), tmp_path / 'config.toml')

    functions = qc_functions()

    assert [qc.name for qc in functions] == ['somalier']
    # The stage's expected outputs are derived from this list, so a restricted run
    # must not claim the outputs of the tools it skipped.
    assert [key for qc in functions for key in qc.outs] == ['somalier']


def test_exome_only_tool_is_selectable_for_exomes(tmp_path: Path):
    set_config(
        config(sequencing_type='exome', tools=['picard_hs_metrics']),
        tmp_path / 'config.toml',
    )

    assert [qc.name for qc in qc_functions()] == ['picard_hs_metrics']


def test_raises_on_a_tool_that_does_not_apply_to_the_sequencing_type(tmp_path: Path):
    # picard_hs_metrics is exome-only, so naming it for a genome is a mistake worth
    # failing on rather than silently running nothing.
    set_config(config(tools=['picard_hs_metrics']), tmp_path / 'config.toml')

    with pytest.raises(ValueError, match='Unknown tool'):
        qc_functions()


def test_raises_on_an_unknown_tool(tmp_path: Path):
    set_config(config(tools=['somalier', 'not_a_tool']), tmp_path / 'config.toml')

    with pytest.raises(ValueError, match='not_a_tool'):
        qc_functions()


def test_skip_qc_still_wins(tmp_path: Path):
    cfg = config(tools=['somalier'])
    cfg.workflow.skip_qc = True
    set_config(cfg, tmp_path / 'config.toml')

    assert qc_functions() == []
