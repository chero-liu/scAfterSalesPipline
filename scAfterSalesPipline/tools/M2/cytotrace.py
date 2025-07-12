def shell_script(self):
    shell_script_contenr = f"""
source /home/liuxuan/miniconda3/bin/activate /gpfs/oe-software/conda_envs/scrna_envs/cytoTRACE

Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M3/cytotrace/CytoTrace_v1.2.R \
    -i data_ob_v3.rds   \
    -o ./ \
    -g clusters \
    -b FALSE

"""
