def find_decoding_error_line(filepath: str, encoding='utf-8'):
    """
    Reads a file and reports the line number of the first UnicodeDecodeError.
    
    This is a diagnostic tool to help pinpoint encoding errors in files.
    """
    try:
        with open(filepath, 'rb') as f: # Open in binary mode
            for i, line_bytes in enumerate(f, 1):
                try:
                    line_bytes.decode(encoding)
                except UnicodeDecodeError as e:
                    print(f"--- Encoding Error Detected in {filepath} ---")
                    print(f"Error on line: {i}")
                    print(f"Original error: {e}")
                    print(f"Problematic line (raw bytes): {line_bytes!r}")
                    # Try to decode with a forgiving encoding to see the content
                    try:
                        line_content = line_bytes.decode('latin-1')
                        print(f"Line content (decoded as latin-1): '{line_content.strip()}'")
                    except:
                        pass # Should not fail, but just in case
                    print("--- End of Report ---")
                    return i # Return the line number
    except FileNotFoundError:
        print(f"File not found: {filepath}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    return None # No error found

# Add this to the end of energy.py
if __name__ == '__main__':
    print("Running energy check as a standalone script...")
    # Check()  # You can comment out the normal run

    print("\nDiagnosing potential encoding issues in 'output.ok'...")
    find_decoding_error_line('output.ok')
