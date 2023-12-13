#!/usr/bin/env python3

import logging
from .cli import mitopy
from ._version import __version__


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        handlers=[logging.StreamHandler()],
    )

    mitopy()


if __name__ == "__main__":
    main()
