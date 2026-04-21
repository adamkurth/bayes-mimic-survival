#!/bin/bash
# Quick check: largest 30 files in repo (excluding .git)
du -sh $(find . -type f ! -path './.git/*' 2>/dev/null) 2>/dev/null | sort -rh | head -30
