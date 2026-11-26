"""Spherical LTE model atmospheres for cool giant stars."""

__all__ = ["initialize", "main", "restart_temperatures", "update_boundary", "write_checkpoint"]


def __getattr__(name):
    if name in __all__:
        from . import cli

        return getattr(cli, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
