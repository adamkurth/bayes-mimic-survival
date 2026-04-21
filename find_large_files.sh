#!/bin/bash

# Script to find large files in the repository before committing to GitHub
# Usage: ./find_large_files.sh [size_threshold_MB]
# Default threshold: 10 MB

THRESHOLD_MB=${1:-10}
THRESHOLD_BYTES=$((THRESHOLD_MB * 1024 * 1024))

echo "================================"
echo "Finding Large Files in Repository"
echo "================================"
echo "Threshold: ${THRESHOLD_MB} MB"
echo ""

echo "--- Large files in working directory ---"
find . -type f -size +${THRESHOLD_MB}M ! -path './.git/*' | while read file; do
    size_bytes=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null)
    size_mb=$(echo "scale=2; $size_bytes / 1048576" | bc)
    printf "%-60s %8.2f MB\n" "$file" "$size_mb"
done | sort -k2 -rn

echo ""
echo "--- Largest 20 files overall ---"
find . -type f ! -path './.git/*' -exec ls -lh {} \; | awk '{print $5, $9}' | sort -h -r | head -20

echo ""
echo "--- Large files by extension ---"
find . -type f ! -path './.git/*' | while read file; do
    ext="${file##*.}"
    stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null
done | awk '{sum+=$1} END {for (ext in sum) print ext, sum[ext]}' 

echo ""
echo "--- Git object analysis ---"
echo "Largest objects in .git/objects:"
du -sh .git/objects/* 2>/dev/null | sort -rh | head -10

echo ""
echo "--- Recommended .gitignore patterns ---"
echo "Consider adding these to .gitignore:"
find . -type f -size +${THRESHOLD_MB}M ! -path './.git/*' | sed 's|^\./||' | sed 's|.*|# &|'
