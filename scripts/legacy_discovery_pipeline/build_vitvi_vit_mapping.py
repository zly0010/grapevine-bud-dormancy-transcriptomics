import gzip, re, json

# Strategy:
# 1. VIT_ IDs from UniProt gene names -> GO terms
# 2. Vitvi -> VIT_ mapping using GTF coordinates or direct pattern matching
# 3. Build final Vitvi -> GO mapping

# Step 1: Build VIT_ -> GO mapping from UniProt
print("Step 1: Building VIT_ -> GO mapping from UniProt...")
with open(r"data/legacy/data/annotation/grape_go_raw_nobom.json", "r") as f:
    data = json.load(f)

vit_to_go = {}  # VIT_ID -> {BP: [(name, id)], MF: [...], CC: [...]}
vit_to_acc = {}  # VIT_ID -> UniProt accession

for r in data:
    acc = r.get("primaryAccession", "")
    
    # Extract VIT_ IDs from gene names
    vit_ids = []
    for g in r.get("genes", []):
        for ol in g.get("orderedLocusNames", []):
            name = ol.get("value", "")
            if name.startswith("VIT_"):
                vit_ids.append(name)
        for sn in g.get("submittedNames", []):
            name = sn.get("value", "")
            if name.startswith("VIT_"):
                vit_ids.append(name)
        if "geneName" in g:
            name = g["geneName"].get("value", "")
            if name.startswith("VIT_"):
                vit_ids.append(name)
    
    if not vit_ids:
        continue
    
    # Extract GO terms
    bp = []
    mf = []
    cc = []
    for xref in r.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "GO":
            go_id = xref.get("id", "")
            props = {p.get("key"): p.get("value") for p in xref.get("properties", [])}
            go_term = props.get("GoTerm", "")
            
            if go_term.startswith("P:"):
                bp.append((go_term[2:], go_id))
            elif go_term.startswith("F:"):
                mf.append((go_term[2:], go_id))
            elif go_term.startswith("C:"):
                cc.append((go_term[2:], go_id))
    
    for vid in set(vit_ids):
        if vid not in vit_to_go:
            vit_to_go[vid] = {"BP": [], "MF": [], "CC": []}
        vit_to_go[vid]["BP"].extend(bp)
        vit_to_go[vid]["MF"].extend(mf)
        vit_to_go[vid]["CC"].extend(cc)
        vit_to_acc[vid] = acc

print(f"  VIT_ IDs with GO: {len(vit_to_go)}")

# Step 2: Build Vitvi -> VIT_ mapping
# Vitvi{chrom}g{number} -> VIT_{chrom}s{scaffold}g{number}
# We need to find which scaffold each gene is on
# Strategy: Use GTF coordinates to match

print("\nStep 2: Loading GTF coordinates...")
f = gzip.open(r"data/legacy/data/annotation/Vitis_viniferaEnsembl62.gtf.gz", "rt")
vit_coords = {}  # Vitis_id -> (chrom, start, end)
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
            vit_coords[gene] = (chrom, start, end)
f.close()
print(f"  Vitis genes with coords: {len(vit_coords)}")

# Step 3: Build VIT_ -> Vitis mapping using UniProt gene names
# VIT_ IDs have scaffold info, Vitis IDs don't
# We need to match them by genomic position or other means

# Actually, let's try a different approach:
# The Ensembl 62 GTF has Vitis IDs
# The UniProt entries have VIT_ IDs in gene names
# Let's check if any UniProt entries have BOTH VIT_ and Vitis names

print("\nStep 3: Checking for cross-references between VIT_ and Vitis...")
vit_to_vitis = {}  # VIT_ -> Vitis
for r in data:
    vit_ids = []
    vitis_ids = []
    for g in r.get("genes", []):
        for ol in g.get("orderedLocusNames", []):
            name = ol.get("value", "")
            if name.startswith("VIT_"):
                vit_ids.append(name)
            elif name.startswith("Vitis"):
                vitis_ids.append(name)
        for sn in g.get("submittedNames", []):
            name = sn.get("value", "")
            if name.startswith("VIT_"):
                vit_ids.append(name)
            elif name.startswith("Vitis"):
                vitis_ids.append(name)
        if "geneName" in g:
            name = g["geneName"].get("value", "")
            if name.startswith("VIT_"):
                vit_ids.append(name)
            elif name.startswith("Vitis"):
                vitis_ids.append(name)
    
    for vid in set(vit_ids):
        for vsid in set(vitis_ids):
            vit_to_vitis[vid] = vsid

print(f"  VIT_ -> Vitis mappings: {len(vit_to_vitis)}")
if vit_to_vitis:
    sample = list(vit_to_vitis.items())[:5]
    for k, v in sample:
        print(f"    {k} -> {v}")

# Step 4: Build Vitvi -> VIT_ mapping
# Vitvi{chrom}g{number} where chrom is 0-based
# VIT_{chrom}s{scaffold}g{number} where chrom is 1-based
# The gene number in Vitvi might correspond to the gene number in VIT_

print("\nStep 4: Building Vitvi -> VIT_ mapping...")

# Load Vitvi IDs
with open(r"data/legacy/results_corrected/07_wgcna_fixed/GSE124820_vst_fixed.txt", "r") as f:
    vitvi_ids = f.readline().split("\t")[1:]
print(f"  Vitvi IDs: {len(vitvi_ids)}")

# Try direct mapping: Vitvi{chrom}g{number} -> VIT_{chrom+1}s* g{number}
# This requires knowing the scaffold for each gene
# Let's try matching by gene number within each chromosome

# First, group VIT_ IDs by chromosome (1-based)
vit_by_chrom = {}
for vid in vit_to_go:
    m = re.match(r"VIT_(\d+)s(\d+)g(\d+)", vid)
    if m:
        chrom = int(m.group(1))
        scaffold = int(m.group(2))
        gene_num = int(m.group(3))
        if chrom not in vit_by_chrom:
            vit_by_chrom[chrom] = {}
        # Key by gene number (across all scaffolds in this chromosome)
        # This is tricky because gene numbers restart per scaffold
        # Let's just store all VIT_ IDs per chromosome
        if gene_num not in vit_by_chrom[chrom]:
            vit_by_chrom[chrom][gene_num] = []
        vit_by_chrom[chrom][gene_num].append(vid)

print(f"  VIT_ IDs by chromosome: {len(vit_by_chrom)}")
for chrom in sorted(vit_by_chrom.keys()):
    print(f"    Chrom {chrom}: {len(vit_by_chrom[chrom])} gene numbers")
