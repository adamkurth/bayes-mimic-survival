#!/bin/bash

# Script to identify and download missing MIMIC-III files

DOWNLOAD_DIR="./physionet.org/files/mimiciii/1.4"
REMOTE_URL="https://physionet.org/files/mimiciii/1.4/"
INDEX_FILE="$DOWNLOAD_DIR/index.html"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}=== MIMIC-III Download Checker ===${NC}\n"

# List of expected files (extracted from index.html)
EXPECTED_FILES=(
    "ADMISSIONS.csv.gz:2525254"
    "CALLOUT.csv.gz:1205036"
    "CAREGIVERS.csv.gz:49529"
    "CHARTEVENTS.csv.gz:4287004212"
    "CPTEVENTS.csv.gz:4971247"
    "DATETIMEEVENTS.csv.gz:55042057"
    "DIAGNOSES_ICD.csv.gz:4717463"
    "DRGCODES.csv.gz:1750041"
    "D_CPT.csv.gz:3951"
    "D_ICD_DIAGNOSES.csv.gz:284950"
    "D_ICD_PROCEDURES.csv.gz:75848"
    "D_ITEMS.csv.gz:187922"
    "D_LABITEMS.csv.gz:11492"
    "ICUSTAYS.csv.gz:1990193"
    "INPUTEVENTS_CV.csv.gz:422105376"
    "INPUTEVENTS_MV.csv.gz:150909733"
    "LABEVENTS.csv.gz:335843735"
    "MICROBIOLOGYEVENTS.csv.gz:7612452"
    "NOTEEVENTS.csv.gz:1165661335"
    "OUTPUTEVENTS.csv.gz:58436520"
    "PATIENTS.csv.gz:571615"
    "PRESCRIPTIONS.csv.gz:103492087"
    "PROCEDUREEVENTS_MV.csv.gz:7814321"
    "PROCEDURES_ICD.csv.gz:1795004"
    "SERVICES.csv.gz:1156392"
    "TRANSFERS.csv.gz:5479949"
    "LICENSE.txt:2542"
    "README.md:64"
    "SHA256SUMS.txt:2340"
)

# Initialize counters
MISSING_COUNT=0
COMPLETE_COUNT=0
INCOMPLETE_COUNT=0
MISSING_FILES=()
MISSING_SIZES=()

echo -e "${YELLOW}Checking which files are missing or incomplete...${NC}\n"

# Check each expected file
for FILE_INFO in "${EXPECTED_FILES[@]}"; do
    FILENAME="${FILE_INFO%:*}"
    EXPECTED_SIZE="${FILE_INFO#*:}"
    FILEPATH="$DOWNLOAD_DIR/$FILENAME"

    if [ ! -f "$FILEPATH" ]; then
        # File doesn't exist
        echo -e "${RED}✗ MISSING${NC}: $FILENAME ($(numfmt --to=iec-i --suffix=B $EXPECTED_SIZE 2>/dev/null || echo $EXPECTED_SIZE bytes))"
        MISSING_FILES+=("$FILENAME")
        MISSING_SIZES+=("$EXPECTED_SIZE")
        ((MISSING_COUNT++))
    else
        # File exists - check size
        ACTUAL_SIZE=$(stat -f%z "$FILEPATH" 2>/dev/null || stat -c%s "$FILEPATH" 2>/dev/null)
        if [ "$ACTUAL_SIZE" -eq "$EXPECTED_SIZE" ]; then
            echo -e "${GREEN}✓ COMPLETE${NC}: $FILENAME"
            ((COMPLETE_COUNT++))
        else
            # Incomplete download
            echo -e "${YELLOW}⚠ INCOMPLETE${NC}: $FILENAME ($(numfmt --to=iec-i --suffix=B $ACTUAL_SIZE 2>/dev/null || echo $ACTUAL_SIZE bytes) / $(numfmt --to=iec-i --suffix=B $EXPECTED_SIZE 2>/dev/null || echo $EXPECTED_SIZE bytes))"
            MISSING_FILES+=("$FILENAME")
            MISSING_SIZES+=("$EXPECTED_SIZE")
            ((INCOMPLETE_COUNT++))
        fi
    fi
done

echo ""
echo -e "${YELLOW}--- Summary ---${NC}"
echo -e "Complete files:    ${GREEN}$COMPLETE_COUNT${NC}"
echo -e "Incomplete files:  ${YELLOW}$INCOMPLETE_COUNT${NC}"
echo -e "Missing files:     ${RED}$MISSING_COUNT${NC}"
echo -e "Total to download: ${RED}$((MISSING_COUNT + INCOMPLETE_COUNT))${NC}"

# Calculate total missing size
TOTAL_MISSING_SIZE=0
for SIZE in "${MISSING_SIZES[@]}"; do
    TOTAL_MISSING_SIZE=$((TOTAL_MISSING_SIZE + SIZE))
done
echo -e "Total missing size: $(numfmt --to=iec-i --suffix=B $TOTAL_MISSING_SIZE 2>/dev/null || echo $TOTAL_MISSING_SIZE bytes)"

# Ask user if they want to download
if [ $((MISSING_COUNT + INCOMPLETE_COUNT)) -gt 0 ]; then
    echo ""
    read -p "Download missing files? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo -e "\n${YELLOW}Starting download...${NC}\n"

        # Build wget command for missing files
        WGET_URLS=""
        for FILENAME in "${MISSING_FILES[@]}"; do
            WGET_URLS="$WGET_URLS $REMOTE_URL$FILENAME"
        done

        # Run wget with resume capability
        wget -c --tries=3 --waitretry=5 --user adamkurth --ask-password $WGET_URLS -P "$DOWNLOAD_DIR"

        echo -e "\n${GREEN}Download complete!${NC}"
        echo -e "Re-running verification...\n"

        # Re-check files
        exec "$0"
    else
        echo "Download cancelled."
    fi
else
    echo -e "\n${GREEN}All files are complete!${NC}"
fi
