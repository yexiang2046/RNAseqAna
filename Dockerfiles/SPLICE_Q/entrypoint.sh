#!/bin/bash

# Verify ps command is available
ps --version >/dev/null 2>&1

# Execute the command passed to the container
exec "$@" 