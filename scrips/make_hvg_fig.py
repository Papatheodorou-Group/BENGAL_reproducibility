import click
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import os


@click.command()
@click.option('--batch_key', type=str, default=None, help="Batch key in identifying HVG and harmony integration")
@click.option('--dir_now', type=str, default=None, help="Dir name")
@click.option('--task_now', type=str, default=None, help="Metadata name")
@click.option('--name_now', type=str, default=None, help="Name name")

def run_hvg(batch_key, dir_now, task_now, name_now):

    os.chdir(f"/nfs/research/irene/ysong/RESULTS/NEXTFLOW/CROSS-SPECIES_test/{dir_now}/results/h5ad_homology_concat")

    ad_oo = sc.read_h5ad(f"{task_now}_one2one_only.h5ad")
    sc.pp.normalize_total(ad_oo, target_sum=1e4)
    sc.pp.log1p(ad_oo)
    sc.pp.highly_variable_genes(ad_oo, batch_key=batch_key)

    ad_he = sc.read_h5ad(f"{task_now}_many_higher_expr.h5ad")
    sc.pp.normalize_total(ad_he, target_sum=1e4)
    sc.pp.log1p(ad_he)
    sc.pp.highly_variable_genes(ad_he, batch_key=batch_key)

    ad_sh = sc.read_h5ad(f"{task_now}_many_higher_homology_conf.h5ad")
    sc.pp.normalize_total(ad_sh, target_sum=1e4)
    sc.pp.log1p(ad_sh)
    sc.pp.highly_variable_genes(ad_sh, batch_key=batch_key)
    

    set1 = set(ad_oo.var.loc[ad_oo.var.highly_variable_intersection, :].index)
    set1= set([i.split('-', 1)[0]  for i in set1])
    
    set2 = set(ad_he.var.loc[ad_he.var.highly_variable_intersection, :].index)
    set2= set([i.split('-', 1)[0]  for i in set2])
    
    set3 = set(ad_sh.var.loc[ad_sh.var.highly_variable_intersection, :].index)
    set3= set([i.split('-', 1)[0]  for i in set3])

    venn3([set1, set2, set3], ('O2O', 'HE', 'SH'))
    plt.title(f"{name_now} HVGs")
    plt.savefig(f"/nfs/research/irene/ysong/RESULTS/NEXTFLOW/CROSS-SPECIES_test/count_HVGs_overlap/{name_now}.png", bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    run_hvg()