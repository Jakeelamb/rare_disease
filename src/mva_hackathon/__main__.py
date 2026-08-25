"""Run the auditable command suite without requiring an editable installation."""

from __future__ import annotations

from .cli import app

if __name__ == "__main__":
    app()
