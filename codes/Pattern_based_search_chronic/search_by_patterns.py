import json
import re

def extract_patterns_from_file(filename):
    """Extract unique patterns matching immunity/immune/immuno from file"""
    patterns = set()
    
    # Read file line by line
    with open(filename, 'r') as f:
        for line in f:
            # Case-insensitive search for immunity/immune/immuno patterns
            if re.search(r'immunity|immune/|immuno', line, re.IGNORECASE):
                # Clean up the line - remove leading/trailing whitespace
                cleaned_line = line.strip()
                if cleaned_line:
                    patterns.add(cleaned_line)
    
    return sorted(list(patterns))

def search_records_by_patterns(filename, patterns):
    """Search JSON records for those containing any of the patterns"""
    with open(filename, 'r') as f:
        data = json.load(f)
    
    matching_records = []
    
    for record in data:
        # Convert record to string for searching
        record_str = json.dumps(record)
        # Check if any pattern is in this record
        for pattern in patterns:
            if pattern in record_str:
                matching_records.append(record)
                break  # No need to check other patterns for this record
    
    return matching_records

def print_matching_records(records, max_records=5):
    """Print matching records in a readable format"""
    print(f"Found {len(records)} matching records:\n")
    
    for i, record in enumerate(records[:max_records]):
        print(f"=== Record {i+1} ===")
        # Try to extract some key fields for readability
        if 'record' in record:
            rec = record['record']
            accessions = rec.get('accessions', ['N/A'])
            print(f"Accessions: {accessions}")
            
            # Look for title in references
            references = rec.get('references', [])
            for ref in references:
                if 'title' in ref:
                    print(f"Title: {ref['title']}")
                    break
            
            # Look for journal in references
            for ref in references:
                if 'journal' in ref:
                    print(f"Journal: {ref['journal']}")
                    break
        else:
            # Fallback to showing part of the record
            record_preview = json.dumps(record, indent=2)[:500] + "..."
            print(record_preview)
        print()

if __name__ == "__main__":
    filename = "noro_all_info"
    
    # Step 1: Extract patterns
    print("Step 1: Extracting immunity/immune/immuno patterns...")
    patterns = extract_patterns_from_file(filename)
    print(f"Found {len(patterns)} unique patterns:")
    for i, pattern in enumerate(patterns, 1):
        print(f"{i:2d}: {pattern}")
    print()
    
    # Step 2: Search for records containing these patterns
    print("Step 2: Searching for records containing these patterns...")
    matching_records = search_records_by_patterns(filename, patterns)
    print(f"Found {len(matching_records)} records containing the patterns.\n")
    
    # Step 3: Show results
    print_matching_records(matching_records)