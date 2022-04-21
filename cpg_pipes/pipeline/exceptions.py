"""
Pipeline errors.
"""


class PipelineError(Exception):
    """
    Error raised by pipeline stages implementation.
    """


class StageInputNotFound(Exception):
    """
    Thrown when a stage requests input from another stage
    that doesn't exist.
    """
