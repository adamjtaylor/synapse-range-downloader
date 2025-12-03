#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Validate input arguments
# -----------------------------
if [ $# -lt 2 ]; then
  echo "Usage: $0 <ENTITY_ID> <BYTE_RANGE> [OUTPUT_FILE]"
  echo "Example: $0 syn1234567 0-1023 chunk.bin"
  exit 1
fi

ENTITY_ID="$1"
RANGE="$2"
OUTFILE="${3:-range.bin}"  # Default output filename if not provided

# -----------------------------
# Extract Synapse auth token from ~/.synapseConfig
# -----------------------------
# The file contains a line like:
#   authToken=eyJ0eX...
# We cut on '=' and trim whitespace.
TOKEN=$(grep authToken ~/.synapseConfig | cut -d'=' -f2 | xargs)

# Check we actually got a token
if [ -z "$TOKEN" ]; then
  echo "Error: could not read authToken from ~/.synapseConfig"
  exit 1
fi

# -----------------------------
# Step 1: Get entity details to extract file handle ID
# -----------------------------
echo "Fetching entity details for ${ENTITY_ID}..."
ENTITY_JSON=$(curl -s \
  -H "Authorization: Bearer $TOKEN" \
  "https://repo-prod.prod.sagebase.org/repo/v1/entity/${ENTITY_ID}")

# Check the entity type
CONCRETE_TYPE=$(echo "$ENTITY_JSON" | jq -r '.concreteType')
ENTITY_NAME=$(echo "$ENTITY_JSON" | jq -r '.name')

echo "Entity: ${ENTITY_NAME} (${CONCRETE_TYPE})"

# Extract the dataFileHandleId from the entity response
FILE_HANDLE_ID=$(echo "$ENTITY_JSON" | jq -r '.dataFileHandleId')

# Validate that we got a file handle ID
if [ "$FILE_HANDLE_ID" = "null" ] || [ -z "$FILE_HANDLE_ID" ]; then
  echo ""
  echo "Error: Entity does not have a dataFileHandleId"
  echo "This entity is a ${CONCRETE_TYPE}, not a File entity."
  echo ""
  echo "This script requires a File entity (org.sagebionetworks.repo.model.FileEntity)."
  echo "Please provide a Synapse ID for an actual file, not a Project, Folder, or Table."
  exit 1
fi

echo "Found file handle ID: ${FILE_HANDLE_ID}"

# -----------------------------
# Step 2: Request a pre-signed download URL from Synapse REST API
# -----------------------------
# Use the entity file endpoint with redirect=false to get the presigned URL
# This endpoint returns the URL directly as plain text, not JSON
echo "Fetching presigned URL..."
PRESIGNED=$(curl -s \
  -H "Authorization: Bearer $TOKEN" \
  "https://repo-prod.prod.sagebase.org/repo/v1/entity/${ENTITY_ID}/file?redirect=false")

# Validate that we received a URL
if [ -z "$PRESIGNED" ]; then
  echo "Error: could not fetch presigned URL"
  exit 1
fi

echo "Got presigned URL"

# -----------------------------
# Step 3: Download the specified byte range using HTTP Range header
# -----------------------------
# This avoids downloading the entire file and pulls only the slice you want.
echo "Downloading bytes ${RANGE} from entity ${ENTITY_ID} (file handle ${FILE_HANDLE_ID}) → ${OUTFILE}"

curl -L -H "Range: bytes=${RANGE}" "$PRESIGNED" -o "$OUTFILE"

# -----------------------------
# Step 4: Create text versions of the downloaded bytes
# -----------------------------
TXTFILE="${OUTFILE%.bin}.txt"
HEXFILE="${OUTFILE%.bin}.hex"

echo "Creating text versions..."

# Plain text output (useful for FASTQ, text files)
cat "$OUTFILE" > "$TXTFILE"

# Hex dump with ASCII (useful for binary files like TIFF)
xxd "$OUTFILE" > "$HEXFILE"

echo ""
echo "Done."
echo "  Binary file: ${OUTFILE}"
echo "  Text file:   ${TXTFILE} (plain text)"
echo "  Hex file:    ${HEXFILE} (hex dump)"