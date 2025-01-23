#!/usr/bin/env python3
import sys
from anianns.anianns import main
import setproctitle

# Set the process title to a custom name
setproctitle.setproctitle("ANI Ann's")

sys.exit(main())