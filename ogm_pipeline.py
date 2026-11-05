#!/usr/bin/env python3
"""OGM molecule-level repeat distance analysis from Bionano Solve contigs."""
import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

PAIR = re.compile(r'\((\d+)\s*,\s*(\d+)\)')


def table_rows(path):
    path = Path(path)
    fields = None
    with path.open(encoding='utf8', errors='replace') as f:
        for line in f:
            if line.startswith('#h'):
                fields = line[2:].strip().split()
            elif line.startswith('#') or not line.strip():
                continue
            else:
                if fields is None:
                    raise ValueError(f'No #h header in {path}')
                parts = line.rstrip('\r\n').split('\t')
                if len(parts) != len(fields):
                    raise ValueError(f'Wrong number of fields in {path}: {len(parts)} != {len(fields)}')
                yield dict(zip(fields, parts))


def write_tsv(path, columns, records):
    with Path(path).open('w', newline='', encoding='utf8') as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
        w.writeheader()
        w.writerows(records)


def discover_contigs(root):
    root = Path(root).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(root)
    candidates = []
    for folder in [root] + sorted(p for p in root.rglob('contigs') if p.is_dir()):
        if folder.name != 'contigs':
            continue
        if (folder / 'annotation' / 'exp_refineFinal1_merged.xmap').is_file() and (folder / 'exp_refineFinal1' / 'alignmol' / 'merge').is_dir():
            candidates.append(folder)
    candidates = list(dict.fromkeys(candidates))
    if not candidates:
        raise FileNotFoundError(f'No valid contigs directory found beneath {root}. Expected annotation/exp_refineFinal1_merged.xmap and exp_refineFinal1/alignmol/merge/')
    if len(candidates) > 1:
        raise ValueError('Multiple eligible contigs directories found; provide a narrower --input:\n' + '\n'.join(map(str, candidates)))
    return candidates[0]


def load_label_sites(path, contig_ids):
    result = defaultdict(dict)
    for row in table_rows(path):
        cid = int(row['CMapId'])
        if cid in contig_ids and int(row['LabelChannel']) == 1:
            result[cid][int(row['SiteID'])] = float(row['Position'])
    return result


def select_contigs(base, chrom, start, end):
    root = base / 'annotation' / 'exp_refineFinal1_merged'
    ann_xmap = Path(str(root) + '.xmap')
    r_cmap = Path(str(root) + '_r.cmap')
    q_cmap = Path(str(root) + '_q.cmap')
    for item in (ann_xmap, r_cmap, q_cmap):
        if not item.is_file():
            raise FileNotFoundError(item)
    overlapping = defaultdict(list)
    for row in table_rows(ann_xmap):
        if int(row['RefContigID']) != chrom:
            continue
        low, high = sorted((float(row['RefStartPos']), float(row['RefEndPos'])))
        if not (low <= start and high >= end):
            continue
        pairs = [(int(r), int(q)) for r, q in PAIR.findall(row['Alignment'])]
        if pairs:
            overlapping[int(row['QryContigID'])].append((row, pairs))
    if not overlapping:
        raise ValueError(f'No assembly alignment spans chr{chrom}:{start}-{end}')
    ref = load_label_sites(r_cmap, {chrom}).get(chrom, {})
    query = load_label_sites(q_cmap, set(overlapping))
    mol_dir = base / 'exp_refineFinal1' / 'alignmol' / 'merge'
    candidate = {}
    for cid, entries in sorted(overlapping.items()):
        stem = f'exp_refineFinal1_contig{cid}'
        files = {ext: mol_dir / (stem + ext) for ext in ('.xmap', '_q.cmap', '_r.cmap')}
        missing = [p for p in files.values() if not p.is_file()]
        if missing:
            raise FileNotFoundError(f'Molecule alignment files missing for contig {cid}: {missing}')
        mol_ref = load_label_sites(files['_r.cmap'], {cid}).get(cid, {})
        if not mol_ref:
            raise ValueError(f'No labels in molecule alignment reference CMAP of contig {cid}')
        mapped_entries = []
        for row, pairs in entries:
            matches = defaultdict(set)
            for rid, qid in pairs:  # XMAP: (Reference SiteID, Query SiteID)
                if rid in ref and qid in query.get(cid, {}) and qid in mol_ref:
                    matches[rid].add(qid)
            unique = {rid: next(iter(ids)) for rid, ids in matches.items() if len(ids) == 1}
            mapped_entries.append((row, unique))
        candidate[cid] = {'files': files, 'mol_ref': mol_ref, 'entries': mapped_entries}
    return ref, candidate


def resolve_locus(ref, candidates, start, end):
    left = sorted((rid for rid, pos in ref.items() if pos < start), key=lambda rid: -ref[rid])[:100]
    right = sorted((rid for rid, pos in ref.items() if pos > end), key=lambda rid: ref[rid])[:100]
    if not left or not right:
        raise ValueError('No reference flanking sites found in reference CMAP')
    for l in left:
        for r in right:
            selection = {}
            for cid, entry in candidates.items():
                possibilities = [(row, mapping) for row, mapping in entry['entries'] if l in mapping and r in mapping and mapping[l] != mapping[r]]
                if not possibilities:
                    break
                selection[cid] = max(possibilities, key=lambda x: float(x[0]['Confidence']))
            if len(selection) == len(candidates):
                # Nearest pair in first-left then first-right preference; preserve same ref labels across contigs.
                return l, r, selection
    lines = []
    for cid, entry in sorted(candidates.items()):
        mapped = set().union(*(set(mapping) for _, mapping in entry['entries']))
        lines.append(f'contig {cid}: mapped={len(mapped)}, left={len(mapped.intersection(left))}, right={len(mapped.intersection(right))}')
    raise ValueError('No common reference flanking-label pair across candidate contigs; pooling unequal intervals is unsafe.\n' + '\n'.join(lines))


def read_distances_for_contig(entry, cid, left_q, right_q, interval_start, interval_end, out_path):
    """Map both reference labels to one query molecule each; fetch exact CMAP SiteIDs."""
    map_xmap = entry['files']['.xmap']
    molecule_pairs = {}
    ambiguous = set()
    for row in table_rows(map_xmap):
        if int(row['RefContigID']) != cid:
            continue
        lower, upper = sorted((float(row['RefStartPos']), float(row['RefEndPos'])))
        if lower > interval_start or upper < interval_end:
            continue
        mid = int(row['QryContigID'])
        matches = defaultdict(set)
        for rid, qid in ((int(r), int(q)) for r, q in PAIR.findall(row['Alignment'])):
            if rid in (left_q, right_q):
                matches[rid].add(qid)
        if len(matches[left_q]) != 1 or len(matches[right_q]) != 1:
            continue
        a, b = next(iter(matches[left_q])), next(iter(matches[right_q]))
        if a == b:
            continue
        if mid in molecule_pairs and molecule_pairs[mid] != (a, b):
            ambiguous.add(mid)
        else:
            molecule_pairs[mid] = (a, b)
    for mid in ambiguous:
        molecule_pairs.pop(mid, None)
    if not molecule_pairs:
        write_tsv(out_path, ['QryContigID', 'distance'], [])
        return 0, len(ambiguous)
    needed = {mid: set(pair) for mid, pair in molecule_pairs.items()}
    found = defaultdict(dict)
    for row in table_rows(entry['files']['_q.cmap']):
        mid = int(row['CMapId'])
        if mid not in needed or int(row['LabelChannel']) != 1:
            continue
        sid = int(row['SiteID'])
        if sid in needed[mid]:
            found[mid][sid] = float(row['Position'])
    records = []
    for mid, (a, b) in sorted(molecule_pairs.items()):
        if a in found[mid] and b in found[mid]:
            distance = abs(found[mid][a] - found[mid][b])
            if math.isfinite(distance) and distance > 0:
                records.append({'QryContigID': mid, 'distance': int(distance)})  # legacy aln_modify truncation
    write_tsv(out_path, ['QryContigID', 'distance'], records)
    return len(records), len(ambiguous)


def merge_distances(paths, output):
    entries = defaultdict(list)
    for cid, path in paths:
        for row in csv.DictReader(Path(path).open(newline=''), delimiter='\t'):
            entries[int(row['QryContigID'])].append((cid, int(row['distance'])))
    all_rows, dup_rows, merged_rows, conflicts = [], [], [], []
    for mid, hits in sorted(entries.items()):
        unique_values = {d for _, d in hits}
        for cid, distance in hits:
            item = {'QryContigID': mid, 'distance': distance, 'source_contig': cid}
            all_rows.append(item)
            if len(hits) > 1:
                dup_rows.append(item)
        if len(unique_values) != 1:
            conflicts.append((mid, hits))
        else:
            merged_rows.append({'QryContigID': mid, 'distance': hits[0][1], 'source_contigs': ','.join(str(x) for x in sorted({c for c, _ in hits}))})
    write_tsv(output / 'all_distances.tsv', ['QryContigID', 'distance', 'source_contig'], all_rows)
    write_tsv(output / 'duplicated_molecules.tsv', ['QryContigID', 'distance', 'source_contig'], dup_rows)
    if conflicts:
        write_tsv(output / 'conflicting_molecules.tsv', ['QryContigID', 'distances'], [{'QryContigID': mid, 'distances': str(hits)} for mid, hits in conflicts])
        raise ValueError(f'{len(conflicts)} molecule IDs have conflicting distances. See conflicting_molecules.tsv; cannot safely select arbitrary values')
    write_tsv(output / 'merged_distances.tsv', ['QryContigID', 'distance', 'source_contigs'], merged_rows)
    return len(all_rows), len(entries), sum(len(v)>1 for v in entries.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, required=True, help='Bionano raw results dir (search recursively for /contigs)')
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--prefix', default='', help='Optional result label')
    ap.add_argument('--locus', required=True, help='chr4:39343732-39350590 (example)')
    args = ap.parse_args()
    match = re.fullmatch(r'(?:chr)?(\d+):(\d+)-(\d+)', args.locus, flags=re.IGNORECASE)
    if not match:
        ap.error('--locus must be like chr4:3934999-3953000')
    chrom, start, end = map(int, match.groups())
    if start >= end:
        ap.error('Locus start must be < end')
    if args.prefix and not re.fullmatch(r'[A-Za-z0-9._-]+', args.prefix):
        ap.error('--prefix allows only letters, digits, period, underscore, hyphen')
    base = discover_contigs(args.input)
    name = f'chr{chrom}_{start}_{end}'
    if args.prefix:
        name = f'{args.prefix}_{name}'
    out = args.output.expanduser().resolve() / name
    out.mkdir(parents=True, exist_ok=True)
    print(f'[OGM] contigs: {base}', flush=True)
    print(f'[OGM] output: {out}', flush=True)
    ref, candidates = select_contigs(base, chrom, start, end)
    left, right, selections = resolve_locus(ref, candidates, start, end)
    print(f'[OGM] reference labels: {left} ({ref[left]:.1f}), {right} ({ref[right]:.1f}); distance {ref[right]-ref[left]:.1f} bp', flush=True)
    write_tsv(out / 'reference_interval.tsv', ['ref_chrom','ref_start','ref_end','ref_left_label','ref_right_label','ref_left_bp','ref_right_bp','reference_distance_bp','n_contigs'], [{
        'ref_chrom': chrom, 'ref_start': start, 'ref_end': end, 'ref_left_label': left, 'ref_right_label': right, 'ref_left_bp': ref[left], 'ref_right_bp': ref[right], 'reference_distance_bp': ref[right]-ref[left], 'n_contigs': len(selections)
    }])
    manifest = []
    distance_paths = []
    for cid, (row, mapping) in sorted(selections.items()):
        ql, qr = mapping[left], mapping[right]
        pos_l, pos_r = candidates[cid]['mol_ref'][ql], candidates[cid]['mol_ref'][qr]
        low, high = min(pos_l, pos_r), max(pos_l, pos_r)
        tsv = out / f'contig{cid}.distances.tsv'
        n, n_amb = read_distances_for_contig(candidates[cid], cid, ql, qr, math.floor(low), math.ceil(high), tsv)
        print(f'[OGM] contig{cid}: {n} molecule measurements, {n_amb} ambiguous molecule alignments', flush=True)
        distance_paths.append((cid, tsv))
        manifest.append({'contig':cid,'orientation':row['Orientation'],'ref_left_label':left,'ref_right_label':right,'assembly_left_label':ql,'assembly_right_label':qr,'assembly_left_bp':pos_l,'assembly_right_bp':pos_r,'reference_distance_bp':ref[right]-ref[left],'n_measured':n,'distance_tsv':str(tsv)})
    write_tsv(out / 'manifest.tsv', ['contig','orientation','ref_left_label','ref_right_label','assembly_left_label','assembly_right_label','assembly_left_bp','assembly_right_bp','reference_distance_bp','n_measured','distance_tsv'], manifest)
    total, unique, duplicated = merge_distances(distance_paths, out)
    print(f'[OGM] total={total}; unique={unique}; duplicated IDs={duplicated}', flush=True)
    if unique < 3:
        raise ValueError('Fewer than 3 unique molecules; GMM cannot be reliably fit')
    return out


if __name__ == '__main__':
    main()
