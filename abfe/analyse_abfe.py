import os
import shutil
import time
import pandas as pd
from loguru import logger
from alchemlyb.workflows import ABFE
from typing import List, Dict

logger.add("abfe_analysis.log", rotation="10 MB")

class ABFEAnalyzer:
    def __init__(
        self,
        project_root: str,
        abfe_subdirs: List[str],
        ignore_folders: List[str] = None,
        legs_config: Dict[str, int] = None,
        temperature: float = 298,
    ):
        self.project_root = project_root
        self.abfe_subdirs = abfe_subdirs
        self.ignore_folders = ignore_folders or []
        self.legs_config = legs_config or {
            "complex": {"rest": 11, "coul": 10, "vdw": 21},
            "ligand": {"coul": 11, "vdw": 21},
        }
        self.temperature = temperature
        logger.info(f"Initialized ABFEAnalyzer for root: {self.project_root}")

    @staticmethod
    def check_folder(path: str):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Folder not found: {path}")

    @staticmethod
    def create_folder(path: str):
        os.makedirs(path, exist_ok=True)
        logger.info(f"Ensured folder exists: {path}")

    @staticmethod
    def copy_dhdl_files(source_base: str, dest_dir: str, legs: Dict[str, int], molecule: str):
        ABFEAnalyzer.create_folder(dest_dir)
        for leg, num_windows in legs.items():
            for i in range(num_windows):
                idx = str(i).zfill(2)
                src = os.path.join(source_base, f"{leg}.{idx}/prod/dhdl.xvg")
                dst = os.path.join(dest_dir, f"dhdl.{leg}.{idx}.xvg")
                if os.path.exists(src):
                    shutil.copy(src, dst)
                    logger.info(f"Copied {molecule} {leg} {idx}: {src} -> {dst}")
                else:
                    logger.warning(f"Missing file: {src}")

    def analyse_dhdl(self, path: str):
        cwd = os.getcwd()
        os.chdir(path)
        logger.info(f"Running alchemlyb workflow in {path}")

        workflow = ABFE(
            units="kcal/mol",
            software="GROMACS",
            dir="./",
            prefix="dhdl",
            suffix="xvg",
            T=self.temperature,
            outdirectory="./",
        )
        workflow.run(
            skiptime=10,
            uncorr="dhdl",
            threshold=50,
            methods=("MBAR", "BAR", "TI"),
            overlap="O_MBAR.pdf",
            breakdown=True,
            n_bootstraps=50,
            forwrev=10,
            n_jobs=1,
        )
        workflow.summary.to_csv("ABFE_summary.csv")
        workflow.convergence.to_csv("ABFE_convergence.csv")
        logger.info(f"Workflow complete for {path}")
        os.chdir(cwd)

    @staticmethod
    def get_rest_dG(filepath: str) -> float:
        with open(filepath, "r") as f:
            for line in f:
                if line.startswith("dG_off"):
                    return float(line.split(",")[1].split(":")[1])
        raise ValueError("dG_off not found")

    @staticmethod
    def extract_dG_from_summary(filepath: str, leg: str, stage: str) -> List[float]:
        df = pd.read_csv(filepath)
        df.columns = ["PROCESS", "ID", "MBAR", "MBAR_Error", "BAR", "BAR_Error", "TI", "TI_Error"]
        row = df.loc[df["ID"] == leg].squeeze()
        return [f"{stage}_{leg}"] + row.iloc[2:].astype(float).tolist()

    def compile_results(self, analysis_path: str):
        logger.info(f"Compiling results in {analysis_path}")
        os.chdir(analysis_path)
        os.chdir("analysis")
        rest_dG = self.get_rest_dG("../boresch/resgen_output.txt")

        records = [
            self.extract_dG_from_summary("ligand/ABFE_summary.csv", "coul", "ligand"),
            self.extract_dG_from_summary("ligand/ABFE_summary.csv", "vdw", "ligand"),
            ["anal_rest", rest_dG, 0, rest_dG, 0, rest_dG, 0],
            self.extract_dG_from_summary("complex/ABFE_summary.csv", "bonded", "complex"),
            self.extract_dG_from_summary("complex/ABFE_summary.csv", "coul", "complex"),
            self.extract_dG_from_summary("complex/ABFE_summary.csv", "vdw", "complex"),
        ]
        df = pd.DataFrame(
            records,
            columns=["state", "MBAR", "MBAR_Error", "BAR", "BAR_Error", "TI", "TI_Error"],
        )

        totals = {
            "state": "total",
            "MBAR": df[df.state.str.startswith("ligand")]["MBAR"].sum() + df[df.state == "anal_rest"]["MBAR"].sum() -
                    df[df.state.str.startswith("complex")]["MBAR"].sum(),
            "MBAR_Error": df["MBAR_Error"].sum(),
            "BAR": df[df.state.str.startswith("ligand")]["BAR"].sum() + df[df.state == "anal_rest"]["BAR"].sum() -
                   df[df.state.str.startswith("complex")]["BAR"].sum(),
            "BAR_Error": df["BAR_Error"].sum(),
            "TI": df[df.state.str.startswith("ligand")]["TI"].sum() + df[df.state == "anal_rest"]["TI"].sum() -
                  df[df.state.str.startswith("complex")]["TI"].sum(),
            "TI_Error": df["TI_Error"].sum(),
        }
        df = pd.concat([df, pd.DataFrame([totals])], ignore_index=True)
        df.to_csv("abfe_results.csv", index=False)
        logger.info(f"Results written to abfe_results.csv")

    def process_complex(self, abfe_path: str):
        logger.info(f"Processing COMPLEX at {abfe_path}")
        self.check_folder(abfe_path)
        os.chdir(abfe_path)
        self.check_folder("complex")
        self.create_folder("analysis")
        self.copy_dhdl_files("complex", "analysis/complex", self.legs_config["complex"], "complex")
        self.analyse_dhdl("analysis/complex")

    def process_ligand(self, abfe_path: str):
        logger.info(f"Processing LIGAND at {abfe_path}")
        self.check_folder(abfe_path)
        os.chdir(abfe_path)
        self.check_folder("ligand")
        self.copy_dhdl_files("ligand", "analysis/ligand", self.legs_config["ligand"], "ligand")
        self.analyse_dhdl("analysis/ligand")

    @staticmethod
    def create_symlink(src_path: str, dst_path: str):
        if not os.path.exists(dst_path):
            os.symlink(src_path, dst_path)
            logger.info(f"Created symlink: {dst_path} -> {src_path}")
        else:
            logger.info(f"Symlink already exists: {dst_path}")

    def run(self, folders: List[str] = None):
        folders = folders or [
            d for d in os.listdir(self.project_root)
            if os.path.isdir(os.path.join(self.project_root, d)) and d not in self.ignore_folders
        ]
        logger.info(f"Found simulation folders: {folders}")

        # Define reference ligand analysis path (first folder + first replicate)
        reference_folder = folders[0]
        reference_subdir = self.abfe_subdirs[0]
        reference_path = os.path.join(self.project_root, reference_folder, reference_subdir)

        # Analyze ligand only ONCE for reference
        self.process_ligand(reference_path)

        # Analyze complex for all folders and subdirs
        for folder in folders:
            for subdir in self.abfe_subdirs:
                replicate_path = os.path.join(self.project_root, folder, subdir)
                self.process_complex(replicate_path)

        # Symlink ligand analyses for other replicates
        ref_ligand_analysis_path = os.path.join(reference_path, "analysis", "ligand")
        for folder in folders:
            for subdir in self.abfe_subdirs:
                replicate_path = os.path.join(self.project_root, folder, subdir)
                ligand_analysis_dst = os.path.join(replicate_path, "analysis", "ligand")
                if replicate_path != reference_path:
                    self.create_symlink(ref_ligand_analysis_path, ligand_analysis_dst)

        # Compile results for all folders/subdirs
        for folder in folders:
            for subdir in self.abfe_subdirs:
                replicate_path = os.path.join(self.project_root, folder, subdir)
                self.compile_results(replicate_path)
