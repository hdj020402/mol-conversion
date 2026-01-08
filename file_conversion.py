import os
import subprocess
from rdkit import Chem
from rdkit import RDLogger
from tqdm import tqdm

RDLogger.DisableLog('rdApp.*')

class FileConversion:
    """
    File conversion methods collection class
    Provides all file-to-file molecular format conversion functions
    """
    
    @staticmethod
    def xyz_to_sdf(xyz_file: str, sdf_file: str):
        """
        Convert XYZ file to SDF file using Open Babel program (file-to-file conversion)
        
        Args:
            xyz_file (str): Input XYZ file path
            sdf_file (str): Output SDF file path
        """
        openbabel_cmd = f'obabel -ixyz {xyz_file} -osdf -O {sdf_file}'
        subprocess.Popen(openbabel_cmd,
                         shell=True,
                         stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL)

    @staticmethod
    def xyz_to_inchi(xyz_file: str) -> str:
        """
        Convert XYZ file to InChI format using Open Babel program (file-to-string conversion)
        
        Args:
            xyz_file (str): Input XYZ file path
            
        Returns:
            str: InChI format string
        """
        if not os.path.exists(xyz_file):
            raise FileNotFoundError(f'File {xyz_file} does not exist!')
        openbabel_cmd = f'obabel -ixyz {xyz_file} -oinchi --readconformer'

        result = subprocess.run(
            openbabel_cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
            )

        return result.stdout.replace('\n', '')

    @staticmethod
    def xyz_to_inchikey(xyz_file: str) -> str:
        """
        Convert XYZ file to InChIKey format using Open Babel program (file-to-string conversion)
        
        Args:
            xyz_file (str): Input XYZ file path
            
        Returns:
            str: InChIKey format string
        """
        openbabel_cmd = f'obabel -ixyz {xyz_file} -oinchikey --readconformer'

        result = subprocess.run(
            openbabel_cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
            )

        return result.stdout.replace('\n', '')

    @staticmethod
    def xyz_to_smiles(xyz_file: str) -> str:
        """
        Convert XYZ file to SMILES format using Open Babel program (file-to-string conversion)
        
        Args:
            xyz_file (str): Input XYZ file path
            
        Returns:
            str: SMILES format string
        """
        openbabel_cmd = f'obabel -ixyz {xyz_file} -osmi'

        result = subprocess.run(
            openbabel_cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
            )

        return result.stdout.split('\t')[0]

    @staticmethod
    def merge_sdf_files(sdf_list: list[str], output_path: str, need_remove: bool = False) -> str:
        """
        Merge multiple SDF files into one file
        
        Args:
            sdf_list (list[str]): List of SDF file paths
            output_path (str): Output file path
            need_remove (bool): Whether to delete source files after merging
            
        Returns:
            str: Output file path
        """
        print(f'total sdf: {len(sdf_list)}')

        mol_list = []
        name_list = []
        for sdf in tqdm(sdf_list):
            mol_name = os.path.basename(sdf).split('.')[0]
            mol = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0]
            mol_list.append(mol)
            name_list.append(mol_name)

        out_sdf_name = output_path
        writer = Chem.SDWriter(out_sdf_name)
        for i, mol in enumerate(mol_list):
            mol.SetProp('_Name', name_list[i])
            writer.write(mol)

        writer.close()

        if need_remove:
            for f in sdf_list:
                os.remove(f)

        return out_sdf_name

    @staticmethod
    def xyz_to_pdb(xyz_file: str, pdb_file: str):
        """
        Convert XYZ file to PDB file
        
        Args:
            xyz_file (str): Input XYZ file path
            pdb_file (str): Output PDB file path
        """
        openbabel_cmd = f'obabel -ixyz {xyz_file} -opdb -O {pdb_file}'
        subprocess.Popen(openbabel_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    @staticmethod
    def xyz_to_mol2(xyz_file: str, mol2_file: str):
        """
        Convert XYZ file to MOL2 file
        
        Args:
            xyz_file (str): Input XYZ file path
            mol2_file (str): Output MOL2 file path
        """
        openbabel_cmd = f'obabel -ixyz {xyz_file} -omol2 -O {mol2_file}'
        subprocess.Popen(openbabel_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    @staticmethod
    def xyz_to_cif(xyz_file: str, cif_file: str):
        """
        Convert XYZ file to CIF file
        
        Args:
            xyz_file (str): Input XYZ file path
            cif_file (str): Output CIF file path
        """
        openbabel_cmd = f'obabel -ixyz {xyz_file} -ocif -O {cif_file}'
        subprocess.Popen(openbabel_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    @staticmethod
    def xyz_to_pdb_string(xyz_file: str) -> str:
        """
        Convert XYZ file to PDB format string
        
        Args:
            xyz_file (str): Input XYZ file path
            
        Returns:
            str: PDB format string
        """
        openbabel_cmd = f'obabel -ixyz {xyz_file} -opdb'
        result = subprocess.run(openbabel_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.stdout

    @staticmethod
    def xyz_to_mol2_string(xyz_file: str) -> str:
        """
        Convert XYZ file to MOL2 format string
        
        Args:
            xyz_file (str): Input XYZ file path
            
        Returns:
            str: MOL2 format string
        """
        openbabel_cmd = f'obabel -ixyz {xyz_file} -omol2'
        result = subprocess.run(openbabel_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.stdout

    @staticmethod
    def xyz_to_cif_string(xyz_file: str) -> str:
        """
        Convert XYZ file to CIF format string
        
        Args:
            xyz_file (str): Input XYZ file path
            
        Returns:
            str: CIF format string
        """
        openbabel_cmd = f'obabel -ixyz {xyz_file} -ocif'
        result = subprocess.run(openbabel_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.stdout