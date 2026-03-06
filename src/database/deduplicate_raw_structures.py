"""

"""
import os
import argparse

def deduplicate_raw_structures(directory : str):
    """
    by_name
    """
    file_names = set()
    remove_counter = 0
    for file in os.listdir(directory):
        if file.split(".")[0] in file_names:
            # remove file
            os.remove(os.path.join(directory, file))
            remove_counter += 1
        else:
            file_names.add(file.split(".")[0])

    print(f"Removed {remove_counter} structures. Leaving with {len(file_names)} structures.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="Path to the raw structure files.")
    args = parser.parse_args()

    deduplicate_raw_structures(args.directory)