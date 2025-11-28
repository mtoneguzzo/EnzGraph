import pandas as pd
import argparse
from collections import Counter
import csv

def detect_delimiter(file_path):
    """Detect delimiter using csv.Sniffer."""
    with open(file_path, 'r', encoding='utf-8') as f:
        sample = f.read(2048)
        sniffer = csv.Sniffer()
        try:
            dialect = sniffer.sniff(sample)
            return dialect.delimiter
        except csv.Error:
            return ','  # Default to comma

def load_csv_with_fallback(path):
    """Load CSV with multiple fallback strategies."""
    try:
        df = pd.read_csv(path, engine='python')
        if len(df.columns) == 1:  # Fallback if header wasn't split
            df = pd.read_csv(path, sep=",", engine='python')
        if len(df.columns) == 1:  # Still broken? Force header=None
            df = pd.read_csv(path, sep=",", header=None, engine='python')
            print(f"⚠ Warning: Header missing or malformed in {path}. Columns auto-generated.")
    except Exception as e:
        raise ValueError(f"Failed to read {path}: {e}")
    return df

def main(input_csv, context_csv, output_csv, top_percentage):
    # Load files with fallback logic
    df = load_csv_with_fallback(input_csv)
    context_df = load_csv_with_fallback(context_csv)

    # Clean headers if present
    df.columns = df.columns.str.strip()
    context_df.columns = context_df.columns.str.strip()

    # Validate columns
    if "hit_score" not in df.columns or "assemblies" not in df.columns:
        raise ValueError(f"Input CSV missing required columns. Found: {df.columns.tolist()}")
    if "Assembly" not in context_df.columns or "Biome_ID" not in context_df.columns:
        raise ValueError(f"Context CSV missing required columns. Found: {context_df.columns.tolist()}")

    total_rows_original = len(df)

    # Map assemblies to biomes
    assembly_to_biome = dict(zip(context_df["Assembly"], context_df["Biome_ID"]))

    def map_biomes(assemblies_str):
        assemblies = [a.split("|")[0].strip() for a in str(assemblies_str).split(";") if a.strip()]
        biomes = [assembly_to_biome[a] for a in assemblies if a in assembly_to_biome]
        return biomes

    df["biomes"] = df["assemblies"].apply(map_biomes)

    # Filter empty biomes
    df = df[df["biomes"].apply(lambda x: len(x) > 0)]
    rows_after_filtering = len(df)

    # Sort by hit_score
    df_sorted = df.sort_values(by="hit_score", ascending=False)

    # Top percentage
    top_n = int(len(df_sorted) * (top_percentage / 100))
    top_hits = df_sorted.head(top_n)

    # Biomes of interest
    biomes_of_interest = [
        "root: Environmental:Aquatic:Freshwater",
        "root: Environmental: Terrestrial: Soil",
        "root: Environmental:Aquatic:Marine",
        "root: Engineered: Solid waste: Composting",
        "root: Engineered: Lab enrichment:Defined media:Anaerobic media",
        "root: Engineered: Wastewater:Nutrient removal:Dissolved organics (anaerobic)"
    ]
    normalized_interest = [b.replace(" ", "") for b in biomes_of_interest]

    biome_hits = df_sorted[df_sorted["biomes"].apply(
        lambda b_list: isinstance(b_list, list) and any(
            isinstance(b, str) and b.replace(" ", "").startswith(b_interest)
            for b_interest in normalized_interest for b in b_list
        )
    )]

    # Combine subsets
    subset = pd.concat([top_hits, biome_hits])

    # Convert biomes list to string BEFORE deduplication
    subset["biomes"] = subset["biomes"].apply(lambda x: "; ".join(x) if isinstance(x, list) else "")
    subset = subset.drop_duplicates()

    # Save output
    subset.to_csv(output_csv, index=False)

    # Summary report
    unique_biomes = set()
    biome_counts = Counter()
    for b_str in subset["biomes"]:
        for b in b_str.split("; "):
            if b.strip():
                unique_biomes.add(b.strip())
                biome_counts[b.strip()] += 1

    print("\n=== SUMMARY REPORT ===")
    print(f"Total rows in original file: {total_rows_original}")
    print(f"Rows after filtering empty biomes: {rows_after_filtering}")
    print(f"Top {top_percentage}% rows: {len(top_hits)}")
    print(f"Biome matches: {len(biome_hits)}")
    print(f"Unique biomes: {len(unique_biomes)}")
    print("\nTop 10 biomes:")
    for biome, count in biome_counts.most_common(10):
        print(f"{biome}: {count}")
    print(f"\nSubset saved to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subset hits based on top percentage and biomes of interest.")
    parser.add_argument("--input", required=True, help="Path to the main hits CSV file (filtered.csv).")
    parser.add_argument("--context", required=True, help="Path to the assembly-biome mapping CSV file (filtered.context.csv).")
    parser.add_argument("--output", required=True, help="Path to save the subset CSV file.")
    parser.add_argument("--top", type=float, default=10, help="Top percentage of hits to include (default: 10).")

    args = parser.parse_args()
    main(args.input, args.context, args.output, args.top)
