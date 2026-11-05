#!/usr/bin/env python3
"""Render count/density OGM distributions using the original GMM function."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from extract_gaussian_alleles import extract_gaussian_alleles


def plot(values, keep, means, sds, weights, reference, title, out, density):
    kept = values[keep]
    edges = np.histogram_bin_edges(kept, bins=30)
    width = float(edges[1] - edges[0])
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.hist(kept, bins=edges, density=density, alpha=0.4, label=f'Included (n={len(kept)})')
    margin = max(50, np.ptp(values) * .05)
    x = np.linspace(min(values)-margin, max(values)+margin, 1000)
    scale = 1 if density else len(kept)*width
    total = np.zeros_like(x)
    for i, (mean, sd, weight) in enumerate(zip(means, sds, weights), 1):
        component = weight * np.exp(-.5*((x-mean)/sd)**2)/(sd*np.sqrt(2*np.pi))
        total += component
        ax.plot(x, component*scale, '-', lw=1.7, label=f'Component {i}: {mean:.1f} bp')
    #ax.plot(x, total*scale, color='black', lw=1.8, label='Gaussian mixture')
    ax.axvline(reference, color='gray', ls=':', label=f'Reference: {reference:.1f} bp')
    if (~keep).any():
        ax.scatter(values[~keep], np.zeros((~keep).sum()), marker='x', color='red', label=f'Excluded (n={(~keep).sum()})')
    ax.set_xlabel(title)
    ax.set_ylabel('Probability density (bp$^{-1}$)' if density else 'Molecule count')
    ax.set_ylim(bottom=0)
    ax.grid(alpha=.15)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def run_gmm(input_file, reference_file, output_dir):
    df = pd.read_csv(input_file, sep='\t')
    if not {'QryContigID', 'distance'}.issubset(df.columns):
        raise ValueError('Missing QryContigID or distance column')
    if df.QryContigID.duplicated().any():
        raise ValueError('Duplicate molecule IDs in merged distances')
    distances = pd.to_numeric(df['distance'], errors='raise').to_numpy(dtype=float)
    if len(distances) < 3 or not np.isfinite(distances).all():
        raise ValueError('Need at least 3 finite molecule distances')
    reference = pd.read_csv(reference_file, sep='\t')
    if len(reference) != 1 or 'reference_distance_bp' not in reference.columns:
        raise ValueError('reference_interval.tsv must have one row with reference_distance_bp')
    ref = float(reference['reference_distance_bp'].iloc[0])
    if not np.isfinite(ref) or ref <= 0:
        raise ValueError('Invalid reference distance')
    _legacy, model, keep = extract_gaussian_alleles(distances)
    keep = np.asarray(keep, dtype=bool)
    if keep.shape != distances.shape or keep.sum() < 2:
        raise ValueError('Too few retained molecules after DBSCAN')
    order = np.argsort(model.means_.ravel())
    means = model.means_.ravel()[order]
    sds = np.sqrt(model.covariances_.reshape(-1))[order]
    weights = model.weights_.ravel()[order]
    if (sds <= 0).any():
        raise ValueError('Invalid Gaussian variance')
    output_dir.mkdir(parents=True, exist_ok=True)
    alleles = pd.DataFrame({'component':np.arange(1,len(means)+1), 'mean_distance_bp':means,
        'sd_bp':sds, 'weight':weights, 'reference_distance_bp':ref, 'expansion_bp':means-ref})
    alleles.to_csv(output_dir/'allele_summary.tsv', sep='\t', index=False)
    results = df.copy()
    results['GMM_included'] = keep
    results['expansion_bp'] = distances-ref
    results['component'] = pd.Series(pd.NA, index=results.index, dtype='Int64')
    results['posterior_probability'] = np.nan
    probs = model.predict_proba(distances[keep].reshape(-1,1))[:,order]
    results.loc[keep,'component'] = probs.argmax(axis=1)+1
    results.loc[keep,'posterior_probability'] = probs.max(axis=1)
    results.to_csv(output_dir/'molecule_assignment.tsv', sep='\t', index=False)
    retained = distances[keep].reshape(-1,1)
    pd.DataFrame([{'n_input':len(distances), 'n_retained':int(keep.sum()), 'n_outliers':int((~keep).sum()),
      'n_components':len(means), 'reference_distance_bp':ref, 'BIC':model.bic(retained),
      'AIC':model.aic(retained)}]).to_csv(output_dir/'model_summary.tsv', sep='\t', index=False)
    for name, values, mu, r, title in (
        ('distance', distances, means, ref, 'Inter-label distance (bp)'),
        ('expansion', distances-ref, means-ref, 0, 'Reference-relative expansion (bp)'),
    ):
        for mode, density in (('density',True),('count',False)):
            plot(values, keep, mu, sds, weights, r, title, output_dir/f'{name}_distribution_{mode}.png', density)
    print(f'[GMM] {len(distances)} unique molecules; {keep.sum()} included; {(~keep).sum()} excluded; {len(means)} Gaussian components', flush=True)
    print(alleles.to_string(index=False), flush=True)
    print(f'[GMM] figures and TSV: {output_dir}', flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('input', type=Path)
    p.add_argument('--reference', type=Path)
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    run_gmm(a.input, a.reference or a.input.parent/'reference_interval.tsv', a.output or a.input.parent/'gmm')
