import csv
import tkinter as tk
from tkinter import filedialog

def build_primer_dict():
    # Hide the Tkinter root window
    root = tk.Tk()
    root.withdraw()

    # Ask for file
    file_path = filedialog.askopenfilename(
        title="Select Primer CSV File",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
    )

    if not file_path:
        print("No file selected.")
        return {}

    primer_map = {}
    with open(file_path, newline='') as f:
        reader = csv.reader(f)
        header = next(reader)  # Skip header row

        # Extract all unique gene names from the first column
        genes = []
        for row in reader:
            if not row:
                continue
            gene = row[0].strip()
            if gene and gene not in genes:
                genes.append(gene)

    # Generate well names in A1..H12 order
    rows = "ABCDEFGH"
    cols = range(1, 13)
    well_names = [f"{r}{c}" for r in rows for c in cols]

    if len(genes) > len(well_names):
        raise RuntimeError(
            f"Too many primers ({len(genes)}) for available wells ({len(well_names)})."
        )

    # Assign each gene to a well in order
    for gene, well in zip(genes, well_names):
        primer_map[gene] = f'primer_rack.wells_by_name()["{well}"]'

    return primer_map

def build_sample_genes():
    # Hide the Tkinter root window
    root = tk.Tk()
    root.withdraw()

    # Ask for file
    file_path = filedialog.askopenfilename(
        title="Select sample/genes CSV File",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
    )

    if not file_path:
        print("No file selected.")
        return {}
    
    sample_genes = {}
    with open(file_path, newline='') as f:
        reader = csv.reader(f)
        header = next(reader)

        for row in reader:
            sample = row[0].strip()
            genes = [g.strip() for g in row[1:] if g.strip()]
            if genes:  # only add if the list is non-empty
                sample_genes[sample] = genes
            else:
                # optional: log that this sample was skipped
                print(f"Skipping control sample {sample} (no genes)")
    
    return sample_genes


if __name__ == "__main__":
    primer_map = build_primer_dict()
    sample_genes = build_sample_genes()

    if primer_map:
        output_file = "primer_map_output.txt"
        with open(output_file, "w") as f:
            f.write("primer_map = {\n")
            for k, v in primer_map.items():
                f.write(f'    "{k}": {v},\n')
            f.write("}\n")
            f.write("sample_genes = {\n")
            for k, v in sample_genes.items():
                f.write(f'    "{k}": {v},\n')
            f.write("}\n")
            
        print(f"\n✅ Dictionaries written to {output_file}. Open it to copy/paste into your protocol.")

