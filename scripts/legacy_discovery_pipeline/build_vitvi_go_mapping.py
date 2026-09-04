import gzip, re, json

# Strategy: Use Ensembl 62 GTF coordinates + UniProt mapping to build Vitvi->GO
# Step 1: Load GTF coordinates for all Vitis genes
print("Step 1: Loading GTF coordinates...")
f = gzip.open(r"data/legacy/data/annotation/Vitis_viniferaEnsembl62.gtf.gz", "rt")
vit_coords = {}
for line in f:
    if line.startswith("#"):
        continue
    parts = line.split("\t")
    if len(parts) > 8 and parts[2] == "gene":
        m = re.search(r'gene_id "([^"]+)"', parts[8])
        if m:
            gene = m.group(1)
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            vit_coords[gene] = (chrom, start, end, strand)
f.close()
print(f"  Vitis genes with coords: {len(vit_coords)}")

# Step 2: Load UniProt mapping: Vitis -> UniProt
print("Step 2: Loading UniProt mapping...")
f = gzip.open(r"data/legacy/data/annotation/grape_uniprot_mapping.tsv.gz", "rt")
header = f.readline()
vitis_to_uniprot = {}
for line in f:
    parts = line.split("\t")
    if len(parts) >= 4:
        v = parts[0]
        u = parts[3]
        if v.startswith("Vitis"):
            if v not in vitis_to_uniprot:
                vitis_to_uniprot[v] = set()
            vitis_to_uniprot[v].add(u)
f.close()
print(f"  Vitis->UniProt: {len(vitis_to_uniprot)}")

# Step 3: Load UniProt GO: UniProt -> GO terms
print("Step 3: Loading GO annotation...")
go_data = {}  # UniProt -> {BP: [terms], MF: [terms], CC: [terms]}
go_ids = {}  # UniProt -> {BP: [GO:ids], MF: [GO:ids], CC: [GO:ids]}
with open(r"data/legacy/data/annotation/uniprot_grape_go.tsv", "r") as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        entry = parts[0]
        bp_terms = []
        bp_ids = []
        mf_terms = []
        mf_ids = []
        cc_terms = []
        cc_ids = []
        if len(parts) >= 3 and parts[2]:
            for term in parts[2].split("; "):
                term = term.strip()
                if "[GO:" in term:
                    go_id = re.search(r'\[(GO:\d+)\]', term)
                    name = term.split(" [GO:")[0].strip()
                    if go_id:
                        bp_terms.append(name)
                        bp_ids.append(go_id.group(1))
        if len(parts) >= 4 and parts[3]:
            for term in parts[3].split("; "):
                term = term.strip()
                if "[GO:" in term:
                    go_id = re.search(r'\[(GO:\d+)\]', term)
                    name = term.split(" [GO:")[0].strip()
                    if go_id:
                        mf_terms.append(name)
                        mf_ids.append(go_id.group(1))
        if len(parts) >= 5 and parts[4]:
            for term in parts[4].split("; "):
                term = term.strip()
                if "[GO:" in term:
                    go_id = re.search(r'\[(GO:\d+)\]', term)
                    name = term.split(" [GO:")[0].strip()
                    if go_id:
                        cc_terms.append(name)
                        cc_ids.append(go_id.group(1))
        if bp_terms or mf_terms or cc_terms:
            go_data[entry] = {"BP": bp_terms, "MF": mf_terms, "CC": cc_terms}
            go_ids[entry] = {"BP": bp_ids, "MF": mf_ids, "CC": cc_ids}
print(f"  UniProt entries with GO: {len(go_data)}")

# Step 4: Build Vitvi -> GO through Vitvi -> Vitis -> UniProt -> GO
print("Step 4: Building Vitvi -> GO mapping...")

# Load Vitvi IDs from WGCNA
with open(r"data/legacy/results_corrected/07_wgcna_fixed/GSE124820_vst_fixed.txt", "r") as f:
    vitvi_ids = f.readline().split("\t")[1:]
print(f"  Vitvi IDs in WGCNA: {len(vitvi_ids)}")

# Build Vitvi -> Vitis mapping
vitvi_to_vitis = {}
for v in vitvi_ids:
    m = re.match(r"Vitvi(\d+)g(\d+)", v)
    if m:
        chrom = int(m.group(1))
        num = int(m.group(2))
        # Try 0-based same chrom
        vitis0 = f"Vitis{chrom:02d}g{num:05d}"
        if vitis0 in vitis_to_uniprot:
            vitvi_to_vitis[v] = vitis0
            continue
        # Try 1-based chrom+1
        vitis1 = f"Vitis{chrom+1:02d}g{num:05d}"
        if vitis1 in vitis_to_uniprot:
            vitvi_to_vitis[v] = vitis1
            continue

print(f"  Vitvi -> Vitis mapped: {len(vitvi_to_vitis)}")

# Build Vitvi -> GO
vitvi_to_go = {}
for vitvi, vitis in vitvi_to_vitis.items():
    uniprots = vitis_to_uniprot.get(vitis, set())
    for up in uniprots:
        if up in go_data:
            if vitvi not in vitvi_to_go:
                vitvi_to_go[vitvi] = {"BP": set(), "MF": set(), "CC": set(), "BP_ids": set(), "MF_ids": set(), "CC_ids": set()}
            vitvi_to_go[vitvi]["BP"].update(go_data[up]["BP"])
            vitvi_to_go[vitvi]["MF"].update(go_data[up]["MF"])
            vitvi_to_go[vitvi]["CC"].update(go_data[up]["CC"])
            vitvi_to_go[vitvi]["BP_ids"].update(go_ids[up]["BP"])
            vitvi_to_go[vitvi]["MF_ids"].update(go_ids[up]["MF"])
            vitvi_to_go[vitvi]["CC_ids"].update(go_ids[up]["CC"])

print(f"  Vitvi with GO: {len(vitvi_to_go)}")
print(f"  Coverage: {len(vitvi_to_go)}/{len(vitvi_ids)} = {len(vitvi_to_go)/len(vitvi_ids)*100:.1f}%")

# Count GO terms per category
total_bp = sum(len(v["BP"]) for v in vitvi_to_go.values())
total_mf = sum(len(v["MF"]) for v in vitvi_to_go.values())
total_cc = sum(len(v["CC"]) for v in vitvi_to_go.values())
print(f"  Total BP terms: {total_bp}")
print(f"  Total MF terms: {total_mf}")
print(f"  Total CC terms: {total_cc}")
