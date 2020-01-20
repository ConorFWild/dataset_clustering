from typing import NamedTuple, Dict
import logging

logging.basicConfig(filename='luigi.log', level=logging.INFO)

import argparse
from pathlib import Path
import os
import shutil

import pandas as pd

import luigi

luigi.configuration.get_config().set('logging', 'log_level', 'INFO')

import dataset_clustering
# from dataset_clustering.cluster import cluster_datasets_rsa as cluster_datasets
from dataset_clustering.cluster import cluster_datasets_rsa as cluster_datasets
from dataset_clustering.data_handling import copy_datasets
from dataset_clustering.cluster import cluster_datasets_luigi

from mdc3.types.datasets import (parse_pandda_input,
                                 parse_pandda_input_for_regex,
                                 )
from mdc3.functions.utils import res_from_mtz_file

from XChemDB.xchem_db import XChemDB

from dataset_clustering.luigi_lib import (CopyDataset,
                                          AlignDataset,
                                          ClusterDatasets,
                                          NormaliseStructure,
                                          normalise_structure,
                                          )


class ClusterXCDBArgs(NamedTuple):
    root_path: Path
    out_dir: Path
    n_procs: int
    clean_run: bool
    pdb_regex: str
    mtz_regex: str
    structure_factors: str

    # xcdb_path: Path


class ClusterArgs(NamedTuple):
    data_dirs: Path
    out_dir: Path
    processor: str
    n_procs: int


def get_args():
    parser = argparse.ArgumentParser()
    # IO
    parser.add_argument("-r", "--root_path",
                        type=str,
                        help="The directory OF THE ROOT OF THE XCHEM DATABASE",
                        required=True
                        )

    parser.add_argument("-o", "--out_dir",
                        type=str,
                        help="The directory for output and intermediate files to be saved to",
                        required=True
                        )

    parser.add_argument("-n", "--n_procs",
                        type=str,
                        help="Number of processes to start",
                        required=True
                        )

    parser.add_argument("-c", "--clean_run",
                        type=bool,
                        help="Number of processes to start",
                        default=True,
                        )

    parser.add_argument("--mtz_regex",
                        type=str,
                        help="Number of processes to start",
                        default="dimple.mtz",
                        )

    parser.add_argument("--pdb_regex",
                        type=str,
                        help="Number of processes to start",
                        default="dimple.pdb",
                        )
    parser.add_argument("--structure_factors",
                        type=str,
                        help="Number of processes to start",
                        default="FWT,PHWT",
                        )

    return parser


def cat_dfs(dfs: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    df = pd.concat(dfs,
                   ignore_index=True,
                   )
    return df


if __name__ == "__main__":

    cmd_args = get_args().parse_args()

    args = ClusterXCDBArgs(root_path=Path(cmd_args.root_path),
                           out_dir=Path(cmd_args.out_dir),
                           n_procs=int(cmd_args.n_procs),
                           clean_run=bool(cmd_args.clean_run),
                           pdb_regex=cmd_args.pdb_regex,
                           mtz_regex=cmd_args.mtz_regex,
                           structure_factors=cmd_args.structure_factors,
                           )

    cluster_dfs: Dict[str, pd.DataFrame] = {}

    print("System intial models: {}".format(args.root_path))

    system_out_dir = args.out_dir
    system_copied_dir = system_out_dir / "copied"
    system_aligned_dir = system_out_dir / "aligned"
    system_processed_dir = system_out_dir / "processed"

    if args.clean_run:

        try:
            shutil.rmtree(system_out_dir.resolve(),
                          ignore_errors=True,
                          )
        except:
            pass

        os.mkdir(system_out_dir.resolve())
        os.mkdir(system_copied_dir.resolve())
        os.mkdir(system_aligned_dir.resolve())
        os.mkdir(system_processed_dir.resolve())

    print("\tParsing pandda input")
    datasets = parse_pandda_input(Path(args.root_path),
                                  mtz_regex=args.mtz_regex,
                                  pdb_regex=args.pdb_regex
                                  )
    print("\t\tFound {} datasets".format(len(datasets)))

    # print([res_from_mtz_file(d["mtz_path"]).limit for dtag, d in datasets.items()])
    truncated_datasets = {}
    for dtag, d in datasets.items():
        if res_from_mtz_file(d["mtz_path"]).limit > 3.0:
            print("\tRemoving dataset {}: bad resolution".format(dtag))
            continue
        else:
            truncated_datasets[dtag] = d

    print("\tGetting reference")
    reference_dtag = list(truncated_datasets.keys())[0]

    print("\tCopying datasets")
    copy_datasets_tasks = [CopyDataset(dtag=dtag,
                                       dataset_path=truncated_datasets[dtag],
                                       output_path=system_copied_dir / dtag,
                                       )
                           for dtag
                           in truncated_datasets
                           ]

    luigi.build(copy_datasets_tasks,
                workers=20,
                local_scheduler=True,
                )

    print("\tSearching for copied datasets in: {}".format(system_copied_dir))
    copied_datasets = parse_pandda_input(system_copied_dir,
                                         pdb_regex="*.pdb",
                                         mtz_regex="*.mtz",
                                         )
    print("\t\tFound {} datasets after copying".format(len(copied_datasets)))

    # normalise_structure_task = NormaliseStructure(reference_pdb_path=copied_datasets[reference_dtag]["pdb_path"],
    #                                               dtag=reference_dtag,
    #                                               output_path=system_out_dir,
    #                                               )
    #
    # luigi.build([normalise_structure_task])

    normalise_structure(copied_datasets[reference_dtag]["pdb_path"],
                        reference_dtag,
                        output_path=system_out_dir,
                        )

    print("\tAligning initial models")
    # align_dataset_tasks = [AlignDataset(dtag=dtag,
    #                                     reference_dtag=reference_dtag,
    #                                     dataset_path=copied_datasets[dtag],
    #                                     reference_dataset_path=copied_datasets[reference_dtag],
    #                                     output_path=system_aligned_dir / dtag,
    #                                     min_res=3.0,
    #                                     )
    #                        for dtag
    #                        in copied_datasets
    #                        ]
    #
    # luigi.build(align_dataset_tasks,
    #             workers=20,
    #             local_scheduler=True,
    #             )

    align_dataset_tasks = [dataset_clustering.luigi_lib.AlignMapToReference(dtag=dtag,
                                                                            reference_dtag=reference_dtag,
                                                                            dataset_path=copied_datasets[dtag],
                                                                            reference_pdb_path=system_out_dir / "{}_normalised.pdb".format(
                                                                                reference_dtag),
                                                                            output_path=system_aligned_dir / dtag,
                                                                            min_res=3.0,
                                                                            structure_factors=args.structure_factors,
                                                                            )
                           for dtag
                           in copied_datasets
                           ]

    luigi.build(align_dataset_tasks,
                workers=20,
                local_scheduler=True,
                )

    # dtag = list(copied_datasets.keys())[10]
    #
    # dataset_path = copied_datasets[dtag]
    # dataset_clustering.luigi_lib.align_map_to_reference(dtag,
    #                                                     reference_dtag,
    #                                                     dataset_path,
    #                                                     system_out_dir / "{}_normalised.pdb".format(reference_dtag),
    #                                                     system_out_dir,
    #                                                     3.0,
    #                                                     )
    #
    #
    # dtag = "MUREECA-x0540"
    # dataset_path = copied_datasets[dtag]
    # dataset_clustering.luigi_lib.align_map_to_reference(dtag,
    #                                                     reference_dtag,
    #                                                     dataset_path,
    #                                                     system_out_dir / "{}_normalised.pdb".format(reference_dtag),
    #                                                     system_out_dir,
    #                                                     3.0,
    #                                                     )
    #
    # dtag = reference_dtag
    # dataset_path = copied_datasets[dtag]
    # dataset_clustering.luigi_lib.align_map_to_reference(dtag,
    #                                                     reference_dtag,
    #                                                     dataset_path,
    #                                                     system_out_dir / "{}_normalised.pdb".format(reference_dtag),
    #                                                     system_out_dir,
    #                                                     3.0,
    #                                                     )

    # exit()

    print("\tClustering aligned models")
    aligned_ccp4s = parse_pandda_input_for_regex(system_aligned_dir,
                                                 "*.ccp4",
                                                 )

    # clustering_task = ClusterDatasets(ccp4_map_paths=aligned_ccp4s,
    #                                   output_path=system_processed_dir,
    #                                   )
    cluster_df = cluster_datasets_luigi(aligned_ccp4s,
                                        system_processed_dir,
                                        "joblib",
                                        20,
                                        )
    # cluster_df.to_csv(str(output_path / "clustering.csv"))

    # luigi.build([clustering_task])

    cluster_dfs_df = cat_dfs(cluster_dfs)

    cluster_dfs_df.to_csv(str(Path(cmd_args.out_dir) / "clusters.csv"))
