"""Self-contained synthetic fixture test for mapping and deduplication."""
import tempfile
from pathlib import Path
from ogm_pipeline import discover_contigs, select_contigs, resolve_locus, read_distances_for_contig, merge_distances


def write_table(path, cols, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as f:
        f.write('#h\t'+'\t'.join(cols)+'\n')
        for row in rows:
            f.write('\t'.join(map(str,row))+'\n')


def test_synthetic():
    with tempfile.TemporaryDirectory() as d:
        base = Path(d)/'deep'/'Bionano'/'output'/'contigs'
        ann = base/'annotation'/'exp_refineFinal1_merged'
        xcols = ['XmapEntryID','QryContigID','RefContigID','QryStartPos','QryEndPos','RefStartPos','RefEndPos','Orientation','Confidence','HitEnum','QryLen','RefLen','LabelChannel','Alignment']
        cmap_cols = ['CMapId','ContigLength','NumSites','SiteID','LabelChannel','Position']
        write_table(Path(str(ann)+'.xmap'),xcols,[[1,221,4,100,200,1000,2000,'-',20,'',1000,100000,1,'(100,10)(101,9)'],[2,3012,4,100,200,1000,2000,'-',20,'',1000,100000,1,'(100,11)(101,10)']])
        write_table(Path(str(ann)+'_r.cmap'),cmap_cols,[[4,100000,2,100,1,1200],[4,100000,2,101,1,1800]])
        write_table(Path(str(ann)+'_q.cmap'),cmap_cols,[[221,1000,2,9,1,200],[221,1000,2,10,1,800],[3012,1000,2,10,1,200],[3012,1000,2,11,1,800]])
        mol = base/'exp_refineFinal1'/'alignmol'/'merge'
        for cid, ql, qr in [(221,10,9),(3012,11,10)]:
            stem=mol/f'exp_refineFinal1_contig{cid}'
            write_table(Path(str(stem)+'_r.cmap'),cmap_cols,[[cid,1000,2,qr,1,210],[cid,1000,2,ql,1,810]])
            xrows=[]; qrows=[]
            for i in range(20):
                mid = i+100 if cid==221 else i+110
                # IDs are shared for 10 molecules; length consistent across contigs.
                dist=600 + (mid % 5)*20
                xrows.append([i+1, mid,cid,100,200,100,900,'+',20,'',1000,1000,1,f'({ql},1)({qr},2)'])
                qrows.extend([[mid,1200,2,1,1,100],[mid,1200,2,2,1,100+dist]])
            write_table(Path(str(stem)+'.xmap'),xcols,xrows)
            write_table(Path(str(stem)+'_q.cmap'),cmap_cols,qrows)
        assert discover_contigs(Path(d)) == base
        ref, candidates=select_contigs(base,4,1300,1700)
        l,r,sel=resolve_locus(ref,candidates,1300,1700)
        assert (l,r)==(100,101)
        paths=[]
        for cid,(row,mapping) in sorted(sel.items()):
            ql,qr=mapping[l],mapping[r]
            p=Path(d)/f'contig{cid}.distances.tsv'
            n, ambiguous=read_distances_for_contig(candidates[cid],cid,ql,qr,210,810,p)
            assert n==20 and ambiguous==0
            paths.append((cid,p))
        out=Path(d)/'out';out.mkdir()
        total,uniq,dups=merge_distances(paths,out)
        assert (total,uniq,dups)==(40,30,10),(total,uniq,dups)
        assert len((out/'merged_distances.tsv').read_text().splitlines())==31
        print('PASS synthetic test: nested detection, two contigs, label mapping, 40 measurements, 30 unique, 10 duplicates')

if __name__=='__main__': test_synthetic()
