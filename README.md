# Synapse Range Request Downloader

A bash script to download specific byte ranges from Synapse files without downloading the entire file. Useful for examining file headers, checking file types, and sampling large files.

## Features

- Downloads only the bytes you need (saves bandwidth and time)
- Automatically chains Synapse REST API endpoints
- Validates entity types (ensures you're downloading from a File entity)
- Outputs multiple formats: binary, plain text, and hex dump
- Works with any file type in Synapse

## Prerequisites

1. **Synapse authentication**: You must have a valid Synapse account with a configured `~/.synapseConfig` file containing your auth token
2. **Required tools**:
   - `curl` - for HTTP requests
   - `jq` - for JSON parsing
   - `xxd` - for hex dump generation
   - `bash` - to run the script
3. **Optional - For composition detection**:
   - `python3` - for running the composition analysis tool
   - `sourmash` - for k-mer sketching and species detection
   - Install with: `pip install -r requirements.txt`

## Installation

```bash
chmod +x test-syn-get-range.sh
```

## Usage

```bash
./test-syn-get-range.sh <ENTITY_ID> <BYTE_RANGE> [OUTPUT_FILE]
```

### Arguments

- `ENTITY_ID`: Synapse entity ID (e.g., `syn1234567`)
- `BYTE_RANGE`: Byte range to download (e.g., `0-1023` for first 1024 bytes)
- `OUTPUT_FILE`: (Optional) Output filename, defaults to `range.bin`

### Examples

#### Check if a file is TIFF or BigTIFF

Download the first 50 bytes to examine the file header:

```bash
./test-syn-get-range.sh syn50944248 0-49 tiff-check.bin
```

Then check the hex output:
```bash
cat tiff-check.hex
```

- Classic TIFF: `4d4d 002a` (big-endian) or `4949 2a00` (little-endian)
- BigTIFF: `4d4d 002b` (big-endian) or `4949 2b00` (little-endian)

#### Extract and Analyze FASTQ Metadata

Download the first 1MB of a gzipped FASTQ file to analyze without downloading the entire file:

```bash
./test-syn-get-range.sh syn59058208 0-1048575 first_mb.fastq.gz
```

Output:
```
Fetching entity details for syn59058208...
Entity: 508_R1.fastq.gz (org.sagebionetworks.repo.model.FileEntity)
Found file handle ID: 140382289
Fetching presigned URL...
Got presigned URL
Downloading bytes 0-1048575 from entity syn59058208 (file handle 140382289) → first_mb.fastq.gz
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 1024k  100 1024k    0     0  1655k      0 --:--:-- --:--:-- --:--:-- 1654k
Creating text versions...

Done.
  Binary file: first_mb.fastq.gz
  Text file:   first_mb.fastq.gz.txt (plain text)
  Hex file:    first_mb.fastq.gz.hex (hex dump)
```

View the first few FASTQ records:
```bash
gunzip -c first_mb.fastq.gz | head
```

Output:
```
@VH00232:284:AAC2F2YHV:1:1101:48604:1076 1:N:0:GATCAAGGCA+ATTAACAAGG
GGCAGGGCGGCCTCCCGCGCCGCGACGTAGGTGATGATGGACTTCTCGTCGGGCTGAGGGACATCCACGTCTGCAGGGAAGGGCAGCTCAGCGGTGGCTCT
+
CC;CC;;;CCC;-;CCC-CCCC;C-CCCCCCCCCCCCCCC;CC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@VH00232:284:AAC2F2YHV:1:1101:50725:1076 1:N:0:GATCAAGGCA+ATTAACAAGG
GGGCGCCTGGGCCGGCGGCGGGGTCAGACCAATCCCGTTGTTTGATTGCCCAGACTACCCCTGTAGCAGGGTCTCAGTCCCTTTCCTCTGGGGCAGTGGCA
+
CCCCC;;C;C;;;;CCCC;C-CCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCC;CCCC;C
@VH00232:284:AAC2F2YHV:1:1101:54436:1076 1:N:0:GATCAAGGCA+ATTAACAAGG
GTGTGGGTAGGCTGGGGTCGGGCGAGGGACAGGCGGTCATTGAAGAGCTGCCCACAGACGTCACATTTGAAGTACTTCTCATCAGCTTGATTGGCACCACC
```

Extract just the header lines:
```bash
gunzip -c first_mb.fastq.gz | grep "^@" | head
```

Analyze metadata with a Python script (10,000 reads from 1MB sample):
```bash
python fastq_metadata.py first_mb.fastq.gz
```

Output shows comprehensive metadata from just the first 1MB:
```json
{
  "instrument": "VH00232",
  "run_number": "284",
  "flowcell": "AAC2F2YHV",
  "lanes": {
    "1": 10000
  },
  "tiles": {
    "1101": 10000
  },
  "index_pairs": {
    "GATCAAGGCA+ATTAACAAGG": 9749,
    "GATCAAGGCA+ATTACCAAGG": 38,
    "GATCAAGGCA+ATTTACAAGG": 16,
    "GATCAAGGCA+ATTAACTAGG": 64,
    ...
  },
  "read_lengths": {
    "101": 10000
  },
  "gc_fraction_est": 0.5245851485148515,
  "qual_summary": {
    "mean_phred": 33.18361782178218
  },
  "reads_parsed": 10000
}
```

This approach allows you to:
- **Quickly assess file quality** without downloading entire files (which can be 10-100GB)
- **Check sequencing parameters** (instrument, run number, flowcell ID)
- **Validate indices** before full download
- **Estimate GC content and quality scores**

#### Detect Species Composition and Contamination (Advanced)

Use sourmash to identify organism composition and detect cross-contamination in FASTQ files. This uses k-mer sketching and containment metrics to determine what organisms are present.

**Prerequisites**:
```bash
# Install sourmash
pip install -r requirements.txt

# Create references directory (if not exists)
mkdir -p references
```

**Setup Reference Signatures**:

You need reference signatures for the organisms you want to detect. You can either download pre-built signatures or create your own:

```bash
# Option 1: Download pre-built signatures (example for common organisms)
# Download from https://sourmash.readthedocs.io/en/latest/databases.html
# or create your own

# Option 2: Create from genome FASTA files
sourmash sketch dna -p k=31,scaled=1000 human_genome.fa -o references/human.sig
sourmash sketch dna -p k=31,scaled=1000 mouse_genome.fa -o references/mouse.sig
sourmash sketch dna -p k=31,scaled=1000 phix.fa -o references/phix.sig
```

**Usage**:
```bash
# Step 1: Download first 1MB using range request
./test-syn-get-range.sh syn59058208 0-1048575 sample.fastq.gz

# Step 2: Detect composition
python fastq_composition.py sample.fastq.gz
```

**Example Output - Cross-Contamination Detected**:
```json
{
  "source": "sample.fastq.gz",
  "reads_sampled": 10000,
  "composition_estimate": {
    "Human": 0.78,
    "Mouse": 0.15
  },
  "unknown_content": 0.07,
  "is_mixed": true,
  "is_contaminated": true,
  "contamination_warning": "Human/Mouse mixture detected"
}
```

**Interpretation**:
- **Pure Sample**: `{"Human": 0.98}` - Species verified, good quality
- **Cross-Contaminated**: `{"Human": 0.78, "Mouse": 0.15}` - Multiple species detected, flag for review
- **Technical Failure**: `{"Phix": 0.95}` - Sequencer control contamination, run failed
- **Unknown**: `{"unknown_content": 0.90}` - Organism not in reference set, may need additional references

**Use Cases**:
- **PDX Model Verification**: Confirm human tumor cells in mouse host (e.g., 80% Human, 20% Mouse expected)
- **Cell Line Authentication**: Verify species matches expected organism
- **Contamination Screening**: Detect PhiX control contamination or cross-sample contamination
- **Species Verification**: Validate metadata before data submission

**Command Options**:
```bash
# Analyze with custom reference directory
python fastq_composition.py sample.fastq.gz --ref-dir ./my_refs

# Limit number of reads processed
python fastq_composition.py sample.fastq.gz --reads 10000

# Stream directly from URL (bypasses bash script)
python fastq_composition.py --url https://example.com/file.fastq.gz
```

#### Sample a large file

Download the first 10KB to inspect structure:

```bash
./test-syn-get-range.sh syn9876543 0-10239 sample.bin
```

## Output Files

The script generates three output files:

1. **`.bin`** - Raw binary data (the actual bytes downloaded)
2. **`.txt`** - Plain text representation (useful for text files like FASTQ, CSV, etc.)
3. **`.hex`** - Hex dump with ASCII visualization (useful for binary files like TIFF, BAM, etc.)

Example output for `chunk.bin`:
```
chunk.bin  - Raw bytes
chunk.txt  - Plain text
chunk.hex  - Hexadecimal dump
```

## How It Works

The script chains together Synapse REST API endpoints:

### Step 1: Entity ID → File Handle ID
```
GET /repo/v1/entity/{ENTITY_ID}
```
Retrieves entity metadata and extracts `dataFileHandleId`

### Step 2: Entity ID → Presigned URL
```
GET /repo/v1/entity/{ENTITY_ID}/file?redirect=false
```
Gets a temporary presigned S3 URL for downloading

### Step 3: Range Request
```
GET {PRESIGNED_URL}
Header: Range: bytes={RANGE}
```
Downloads only the specified byte range using HTTP Range header

### Step 4: Format Conversion
Generates text and hex versions of the downloaded bytes

## File Type Detection

### TIFF Files
- First 4 bytes indicate byte order and version:
  - `4d4d 002a` = Big-endian Classic TIFF
  - `4949 2a00` = Little-endian Classic TIFF
  - `4d4d 002b` = Big-endian BigTIFF
  - `4949 2b00` = Little-endian BigTIFF

### FASTQ Files
- Text-based format with 4 lines per record:
  ```
  @HEADER_LINE
  SEQUENCE
  +
  QUALITY_SCORES
  ```
- Recommend downloading at least 2000 bytes to capture complete records

### BAM/CRAM Files
- Check magic bytes at the beginning:
  - BAM: `1f8b` (gzip header)
  - CRAM: `4352 414d` ("CRAM" in ASCII)

## Error Handling

The script validates:
- Synapse authentication token exists
- Entity is a File entity (not Project, Folder, or Table)
- Presigned URL was successfully retrieved
- All required tools are available

## Limitations

- Byte ranges may cut off in the middle of records for text files
- Presigned URLs expire after 15 minutes (900 seconds)
- Requires read access to the Synapse entity

## Troubleshooting

### "Error: could not read authToken from ~/.synapseConfig"
- Ensure you're logged in to Synapse
- Check that `~/.synapseConfig` exists and contains a valid token

### "Entity does not have a dataFileHandleId"
- You provided a Project, Folder, or Table ID instead of a File entity
- Navigate to an actual file in Synapse and use its entity ID

### "Error: could not fetch presigned URL"
- Check that you have download permissions for the file
- The file may be in a restricted project

## License

MIT

## Contributing

Feel free to submit issues or pull requests for improvements.
