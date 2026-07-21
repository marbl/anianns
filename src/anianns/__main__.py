#!/usr/bin/env python3
import sys

import setproctitle

from anianns.anianns import main as anianns_main


def main():
    """Console entry point kept import-safe for multiprocessing workers."""
    setproctitle.setproctitle("AniAnns")
    return anianns_main()


if __name__ == "__main__":
    sys.exit(main())
