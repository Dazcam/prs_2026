# Extract European samples from genotypes
import pandas as pd
import sys

def main():
    df = pd.read_excel(snakemake.input[0])

    euro_samples = df[df["ancestry"] == "European"]["sample ID"].astype(str).str.strip()

    if euro_samples.empty:
        print("Error: No European samples found in annotation sheet!", file=sys.stderr)
        sys.exit(1)

    euro_samples.to_csv(snakemake.output[0], index=False, header=False)

    print(f"Extracted {len(euro_samples)} European samples → {snakemake.output[0]}", 
          file=sys.stderr)

if __name__ == "__main__":
    main()
