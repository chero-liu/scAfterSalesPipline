import pandas as pd
import os
import argparse
import sys

parser = argparse.ArgumentParser(description="Group and save data by specified column.")


parser.add_argument("--input", type=str, required=True, help="Input file path.")
parser.add_argument("--outdir", type=str, required=True, help="Output directory path.")
parser.add_argument(
    "--splitby",default="None",
    type=str,
    help="Column name to use for grouping.",
)
parser.add_argument(
    "--sortby",default="None",
    type=str,
    help="Column name to use for sorting. From small to large",
)
parser.add_argument(
    "--top",default="None",
    help="",
)
parser.add_argument(
    "--last",default="None",
    help="",
)
args = parser.parse_args()

def createDir(outdir, rm=False):
    if type(outdir).__name__ == 'list':
        for out in outdir:
            if not os.path.exists(out): 
                os.makedirs(out)
            else:
                if rm:
                    os.system(f"rm -r {out}")
                    os.makedirs(out)
    else:
        if not os.path.exists(outdir): 
            os.makedirs(outdir)
        else:
            if rm:
                os.system(f"rm -r {outdir}")
                os.makedirs(outdir)

createDir(args.outdir, rm=False)
data = pd.read_table(args.input)
if args.splitby != "None":
    data = data.groupby(args.splitby)

    for group_name, group_df in data:
        output_file_path = os.path.join(args.outdir, f"{args.splitby}_{group_name}.xls")
        if args.sortby != "None":
            sorted_df = group_df.sort_values(args.sortby, ascending=True)
        else:
            sorted_df = group_df
        if args.top != "None" and args.last != "None":
            sys.exit("top and last are mutually exclusive")
        if args.top != "None":
            sorted_df = sorted_df.head(int(args.top))
        if args.last != "None":
            sorted_df = sorted_df.tail(int(args.last))
        group_df = sorted_df
        try:
            group_df.to_csv(output_file_path, index=False, sep="\t")
            print(f"Group {group_name} has been saved to {output_file_path}.")
        except Exception as e:
            print(f"Failed to save {output_file_path}: {e}")

else:
    data.to_csv(f"{args.outdir}/all.xls", index=False, sep="\t")
