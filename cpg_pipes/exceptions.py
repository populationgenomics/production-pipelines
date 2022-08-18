"""
Exception classes
"""


class PipelineError(Exception):
    """
    Error raised by pipeline stages implementation.
    """


class StageInputNotFoundError(Exception):
    """
    Thrown when a stage requests input from another stage
    that doesn't exist.
    """


class MetamistError(Exception):
    """
    Error while interacting with Metamist.
    """

    pass
