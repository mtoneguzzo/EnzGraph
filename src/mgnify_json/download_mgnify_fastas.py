import pandas as pd
import requests
from pathlib import Path

# === CONFIGURATION ===
csv_path = Path("/Users/christofdeboom/Library/CloudStorage/OneDrive-KULeuven/IBP/Mora_updatedPipeline/similarity_json/csv_subset_boi/subset_PF10503.csv")  # <-- Update this path
output_dir = Path("downloaded_fastas")
ok_dir = output_dir / "ok"
nf_dir = output_dir / "not_found"
ok_list_file = output_dir / "downloaded.txt"
nf_list_file = output_dir / "not_found.txt"

# === SETUP OUTPUT DIRECTORIES ===
ok_dir.mkdir(parents=True, exist_ok=True)
nf_dir.mkdir(parents=True, exist_ok=True)
ok_list_file.write_text("")
nf_list_file.write_text("")

# === READ MGYP IDs FROM CSV ===
df = pd.read_csv(csv_path)
mgyp_ids = df.iloc[:, 0].dropna().unique()  # Assumes MGYP IDs are in the first column

# === FUNCTION TO DOWNLOAD FASTA ===
def download_fasta(mgyp_id):
    url = f"https://api.esmatlas.com/fetchSequence/{mgyp_id}"
    headers = {'Accept': 'application/json'}
    try:
        response = requests.get(url, headers=headers, timeout=30)
        if response.status_code == 200:
            data = response.json()
            sequence = data.get("sequence", "")
            if sequence:
                fasta_content = f">{mgyp_id}\n" + "\n".join(sequence[i:i+60] for i in range(0, len(sequence), 60))
                with open(ok_dir / f"{mgyp_id}.fasta", "w") as f:
                    f.write(fasta_content)
                with open(ok_list_file, "a") as f:
                    f.write(f"{mgyp_id}\n")
                return True
    except Exception:
        pass
    with open(nf_list_file, "a") as f:
        f.write(f"{mgyp_id}\n")
    return False

# === DOWNLOAD LOOP ===
for i, mgyp_id in enumerate(mgyp_ids, 1):
    print(f"[{i}/{len(mgyp_ids)}] Downloading {mgyp_id}...")
    download_fasta(mgyp_id)