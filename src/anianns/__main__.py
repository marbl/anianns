#!/usr/bin/env python3
from contextlib import redirect_stderr, redirect_stdout
import os
import sys


def _load_entrypoint():
    """Suppress import-time library output when quiet mode was requested."""
    quiet = "-q" in sys.argv[1:] or "--quiet" in sys.argv[1:]
    if quiet:
        with open(os.devnull, "w") as null_stream:
            with redirect_stdout(null_stream), redirect_stderr(null_stream):
                import setproctitle
                from anianns.anianns import main as anianns_main
    else:
        import setproctitle
        from anianns.anianns import main as anianns_main
    return setproctitle, anianns_main


setproctitle, anianns_main = _load_entrypoint()


def main():
    """Console entry point kept import-safe for multiprocessing workers."""
    setproctitle.setproctitle("AniAnns")
    return anianns_main()


if __name__ == "__main__":
    sys.exit(main())
