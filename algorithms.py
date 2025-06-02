import io
from PIL import Image
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys, Draw 
from rdkit.Chem.Draw import SimilarityMaps 
import numpy as np
import pandas as pd
from abc import ABC, abstractmethod
import base64

def check_textarea_input(textarea: str):
    smiles_str = [x.strip() for x in textarea.split(",")]
    for smiles in smiles_str:
        molecule = Chem.MolFromSmiles(smiles, sanitize=False)
        if molecule is None:
            return False
    for element in smiles_str[:]:
        if not element.strip():
            smiles_str.remove(element)
    return smiles_str

class FingerprintGenerator(ABC):
    @abstractmethod
    def get_fingerprint_generator():
        pass

    @abstractmethod
    def generate_fingerprints():
        pass

class RDKitFingerprintGenerator(FingerprintGenerator):
    def get_fingerprint_generator(self, data:dict):
        return AllChem.GetRDKitFPGenerator(minPath=data['min_path'], maxPath=data['max_path'], fpSize=data['fps_rdkit'])
    
    def generate_fingerprints(self, molecules:list, data:dict):
        fpgen = self.get_fingerprint_generator(data)
        fps = [fpgen.GetFingerprint(x) for x in molecules]
        return fps
    
    def generate_fingerprints_with_ao(self, molecules:list, data:dict):
        fpgen = self.get_fingerprint_generator(data)
        fps = []
        additional_outputs = []

        for mol in molecules:
            ao = AllChem.AdditionalOutput()
            ao.CollectBitPaths()
            fp = fpgen.GetFingerprint(mol, additionalOutput=ao)
            fps.append(fp)
            info = ao.GetBitPaths()
            additional_outputs.append(info)
        
        return fps, additional_outputs
    
class AtomPairsFingerprintGenerator(FingerprintGenerator):
    def get_fingerprint_generator(self, data:dict):
        return AllChem.GetAtomPairGenerator(fpSize=data['fps_atompairs'])

    def generate_fingerprints(self, molecules:list, data:dict):
        fpgen = self.get_fingerprint_generator(data)
        fps = [fpgen.GetFingerprint(x) for x in molecules]
        return fps

class MorganFingerprintGenerator(FingerprintGenerator):
    def get_fingerprint_generator(self, data:dict):
        return AllChem.GetMorganGenerator(radius=data['radius'], fpSize=data['fps_morgan'])
    
    def generate_fingerprints(self, molecules:list, data:dict):
        fpgen = self.get_fingerprint_generator(data)
        fps = [fpgen.GetFingerprint(x) for x in molecules]
        return fps
    
    def generate_fingerprints_with_ao(self, molecules:list, data:dict):
        fpgen = self.get_fingerprint_generator(data)
        fps = []
        additional_outputs = []

        for mol in molecules:
            ao = AllChem.AdditionalOutput()
            ao.CollectBitInfoMap()
            fp = fpgen.GetFingerprint(mol, additionalOutput=ao)
            fps.append(fp)
            info = ao.GetBitInfoMap()
            additional_outputs.append(info)
        
        return fps, additional_outputs

class MACCSKeysFingerprintGenerator(FingerprintGenerator):
    def get_fingerprint_generator():
        pass
    
    def generate_fingerprints(self, molecules:list, data:dict):
        fps = [MACCSkeys.GenMACCSKeys(x) for x in molecules]
        return fps
    
class SimilarityStrategy(ABC):
    @abstractmethod
    def generate_similarity_matrix():
        pass

class TanimotoStrategy(SimilarityStrategy):
    def generate_similarity_matrix(self, fingerprints:list):
        matrix = np.zeros((len(fingerprints), len(fingerprints)))
        for i in range(len(fingerprints)):
            matrix[i, :] = DataStructs.BulkTanimotoSimilarity(fingerprints[i], fingerprints)
        return matrix
    
class DiceStrategy(SimilarityStrategy):
    def generate_similarity_matrix(self, fingerprints:list):
        matrix = np.zeros((len(fingerprints), len(fingerprints)))
        for i in range(len(fingerprints)):
            matrix[i, :] = DataStructs.BulkDiceSimilarity(fingerprints[i], fingerprints)
        return matrix

class CosineStrategy(SimilarityStrategy):
    def generate_similarity_matrix(self, fingerprints:list):
        matrix = np.zeros((len(fingerprints), len(fingerprints)))
        for i in range(len(fingerprints)):
            matrix[i, :] = DataStructs.BulkCosineSimilarity(fingerprints[i], fingerprints)
        return matrix
    
class SokalStrategy(SimilarityStrategy):
    def generate_similarity_matrix(self, fingerprints:list):
        matrix = np.zeros((len(fingerprints), len(fingerprints)))
        for i in range(len(fingerprints)):
            matrix[i, :] = DataStructs.BulkSokalSimilarity(fingerprints[i], fingerprints)
        return matrix

class RusselStrategy(SimilarityStrategy):
    def generate_similarity_matrix(self, fingerprints:list):
        matrix = np.zeros((len(fingerprints), len(fingerprints)))
        for i in range(len(fingerprints)):
            matrix[i, :] = DataStructs.BulkRusselSimilarity(fingerprints[i], fingerprints)
        return matrix
    
class KulczynskiStrategy(SimilarityStrategy):
    def generate_similarity_matrix(self, fingerprints:list):
        matrix = np.zeros((len(fingerprints), len(fingerprints)))
        for i in range(len(fingerprints)):
            matrix[i, :] = DataStructs.BulkKulczynskiSimilarity(fingerprints[i], fingerprints)
        return matrix

class McConnaugheyStrategy(SimilarityStrategy):
    def generate_similarity_matrix(self, fingerprints:list):
        matrix = np.zeros((len(fingerprints), len(fingerprints)))
        for i in range(len(fingerprints)):
            matrix[i, :] = DataStructs.BulkMcConnaugheySimilarity(fingerprints[i], fingerprints)
        return matrix
    
class TverskyStrategy(SimilarityStrategy):
    def generate_similarity_matrix(self, fingerprints:list, a: int, b: int):
        matrix = np.zeros((len(fingerprints), len(fingerprints)))
        for i in range(len(fingerprints)):
            matrix[i, :] = DataStructs.BulkTverskySimilarity(fingerprints[i], fingerprints, a, b)
        return matrix
    
class DataFrameGenerator:
    def __init__(self, smiles:list, generation_strategy:str, similarity_strategy:str, data:dict):
        self.smiles = smiles
        self.molecules = [Chem.MolFromSmiles(x) for x in smiles]
        self.generation_strategy = self.get_generation_strategy(generation_strategy)
        self.similarity_strategy = self.get_similarity_strategy(similarity_strategy)
        self.data = data

    def get_generation_strategy(self, strategy_name:str):
        if strategy_name == "RDKit":
            return RDKitFingerprintGenerator()
        elif strategy_name == "Morgan":
            return MorganFingerprintGenerator()
        elif strategy_name == "AtomPairs":
            return AtomPairsFingerprintGenerator()
        elif strategy_name == "MACCS Keys":
            return MACCSKeysFingerprintGenerator()
        else:
            raise ValueError(f"Unknown generation strategy: {strategy_name}")
    
    def get_similarity_strategy(self, strategy_name:str):
        if strategy_name == "Tanimoto":
            return TanimotoStrategy()
        elif strategy_name == "Dice":
            return DiceStrategy()
        elif strategy_name == "Cosine":
            return CosineStrategy()
        elif strategy_name == "Russel":
            return RusselStrategy()
        elif strategy_name == "Sokal":
            return SokalStrategy()
        elif strategy_name == "Kulczynski":
            return KulczynskiStrategy()
        elif strategy_name == "McConnaughey":
            return McConnaugheyStrategy()
        elif strategy_name == "Tversky":
            return TverskyStrategy()
        else:
            raise ValueError(f"Unknown similarity strategy: {strategy_name}") 
    
    def get_data_frame(self):
        fingerprints = self.generation_strategy.generate_fingerprints(self.molecules, self.data)
        if isinstance(self.similarity_strategy, TverskyStrategy):
            similarity_matrix = self.similarity_strategy.generate_similarity_matrix(fingerprints, float(self.data['a']), float(self.data['b']))
        else:
            similarity_matrix = self.similarity_strategy.generate_similarity_matrix(fingerprints)
        df = pd.DataFrame(similarity_matrix, index=self.smiles, columns=self.smiles)
        return df.round(2)
    
    def get_molecule_image(self, smiles:str):
        mol = Chem.MolFromSmiles(smiles)
        canvas = Draw.MolDraw2DCairo(1000, 1000)
        canvas.DrawMolecule(mol)
        canvas.FinishDrawing()
        bio = io.BytesIO(canvas.GetDrawingText())
        image = Image.open(bio)
        return image
    
    def get_fingerprint_indices(self):
        fingerprints = self.generation_strategy.generate_fingerprints(self.molecules, self.data)
        fingerprint_dict = dict()
        for i, fingerprint in enumerate(fingerprints):
            fingerprint_on_bits = list(fingerprint.GetOnBits())
            fingerprint_dict[self.smiles[i]] = fingerprint_on_bits
        return fingerprint_dict
    
    def get_fingerprint_bit_image(self, smiles:str, bit: str):
        mol = Chem.MolFromSmiles(smiles)

        if isinstance(self.generation_strategy, MorganFingerprintGenerator):
            fpgen = self.generation_strategy.get_fingerprint_generator(self.data)
            ao = AllChem.AdditionalOutput()
            ao.CollectBitInfoMap()
            fp = fpgen.GetFingerprint(mol, additionalOutput=ao)
            bi = ao.GetBitInfoMap()
            image = Draw.DrawMorganBit(mol, int(bit), bi, useSVG=True)
            svg_bytes = image.encode('utf-8')
            encoded = base64.b64encode(svg_bytes).decode('utf-8')
            return f"data:image/svg+xml;base64,{encoded}"
        elif isinstance(self.generation_strategy, RDKitFingerprintGenerator):
            fpgen = self.generation_strategy.get_fingerprint_generator(self.data)
            ao = AllChem.AdditionalOutput()
            ao.CollectBitPaths()
            fp = fpgen.GetFingerprint(mol, additionalOutput=ao)
            bi = ao.GetBitPaths()
            image = Draw.DrawRDKitBit(mol, int(bit), bi, useSVG=True)
            svg_bytes = image.encode('utf-8')
            encoded = base64.b64encode(svg_bytes).decode('utf-8')
            return f"data:image/svg+xml;base64,{encoded}"
        else:
            raise ValueError("Fingerprint bit images can be only generated for RDKit and Morgan fingerprints.")
        
    def get_similarity_metric(self):
        if isinstance(self.similarity_strategy, TanimotoStrategy):
            return DataStructs.TanimotoSimilarity
        elif isinstance(self.similarity_strategy, DiceStrategy):
            return DataStructs.DiceSimilarity
        elif isinstance(self.similarity_strategy, CosineStrategy):
            return DataStructs.CosineSimilarity
        elif isinstance(self.similarity_strategy, RusselStrategy):
            return DataStructs.RusselSimilarity
        elif isinstance(self.similarity_strategy, SokalStrategy):
            return DataStructs.SokalSimilarity
        elif isinstance(self.similarity_strategy, KulczynskiStrategy):
            return DataStructs.KulczynskiSimilarity
        elif isinstance(self.similarity_strategy, McConnaugheyStrategy):
            return DataStructs.McConnaugheySimilarity
        
    def get_similarity_map(self, smiles1:str, smiles2:str):
        canvas = Draw.MolDraw2DCairo(800, 550)
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2) 
        if isinstance(self.generation_strategy, AtomPairsFingerprintGenerator):
            _, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol1, mol2, lambda m, idx: SimilarityMaps.GetAPFingerprint(m, atomId=idx, fpType='bv'), canvas, metric=self.get_similarity_metric())
            canvas.FinishDrawing()
            bio = io.BytesIO(canvas.GetDrawingText())
            image = Image.open(bio)
            return image
        elif isinstance(self.generation_strategy, MorganFingerprintGenerator):
            _, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol1, mol2, lambda m, idx: SimilarityMaps.GetMorganFingerprint(m, atomId=idx, fpType='bv', radius=int(self.data['radius'])), canvas, metric=self.get_similarity_metric())
            canvas.FinishDrawing()
            bio = io.BytesIO(canvas.GetDrawingText())
            image = Image.open(bio)
            return image