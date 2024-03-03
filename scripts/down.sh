#!/bin/bash

# Base URL for accessing ENA API
ENA_BASE_URL="https://www.ebi.ac.uk/ena/portal/api"
# File containing accessions to process
ACCESSION_FILE="accessions.txt"
# Directory to save downloaded files
OUTPUT_DIR="downloaded_files"

mkdir -p "$OUTPUT_DIR"

while IFS= read -r accession; do
    DOWNLOAD_URL="$ENA_BASE_URL/filereport?accession=${accession}&result=read_run&fields=fastq_ftp"
    # Use awk to extract the FTP URLs, split by semicolon, and remove leading/trailing whitespaces
    FTP_URLS=$(curl -s "$DOWNLOAD_URL" | awk -F'\t' 'NR>1 {gsub(/; /, ";"); print $1}' | tr -d ' ')
    # Print the download url
    echo "FTP_URLS: $FTP_URLS"
    # Loop through each FTP URL and download the file
    IFS=';' read -ra FTP_URL_ARRAY <<< "$FTP_URLS"
    for FTP_URL in "${FTP_URL_ARRAY[@]}"; do
        FILENAME=$(basename "$FTP_URL")
        # Download file with curl, retrying if necessary
        curl -O -J -L -# --retry 3 --retry-delay 5 --retry-max-time 60 "$FTP_URL" -o "$OUTPUT_DIR/$FILENAME"
    done

done < "$ACCESSION_FILE"
