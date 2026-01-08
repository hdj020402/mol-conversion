# 分子转换工具 (Molecular Conversion Tools)

一个专业的分子结构格式转换工具集，支持XYZ格式与其他多种常见化学格式之间的转换。

## 功能特性

- **双模式转换**: 支持文件到文件转换和内存到内存转换
- **多种格式支持**: XYZ ↔ SDF, MOL, InChI, InChIKey, SMILES, PDB, MOL2, CIF
- **GCN编码**: 提供基于图卷积网络的分子编码功能
- **分子标准化**: 自动化的分子坐标标准化处理
- **批处理支持**: 高效的批量分子转换能力
- **面向对象设计**: 所有功能都集成在类中，提供清晰的API接口

## 项目结构

```
mol_conversion/
├── __init__.py               # 包初始化文件
├── setup.py                  # 安装脚本（向后兼容）
├── pyproject.toml           # 现代包配置
├── mol_conversion.py        # 文件转换方法（subprocess调用Open Babel）
├── memory_conversion.py     # 内存转换方法（使用Open Babel Python API）
├── utils.py                 # 分子分析工具（GCN编码、标准化等）
├── mol_conversion.ipynb     # 使用示例Jupyter笔记本
├── IFLOW.md                 # 项目上下文文档
└── README.md                # 本文件
```

## 安装

### 方法1：开发模式安装（推荐）

```bash
# 进入项目目录
cd mol_conversion

# 开发模式安装（可编辑模式）
pip install -e .

# 或者安装完整依赖
pip install -e ".[dev]"
```

### 方法2：直接安装

```bash
# 如果发布到PyPI后
pip install mol-conversion
```

### 方法3：Git仓库安装

```bash
pip install git+https://github.com/your-username/mol-conversion.git
```

### 依赖要求

**系统依赖（重要）:**
- **Open Babel**: 必须通过conda安装，不能通过pip安装
- Python 3.8+

**Python依赖:** (自动安装)
- numpy >= 1.20.0
- rdkit >= 2020.03.5
- tqdm >= 4.60.0

**重要：Open Babel安装（必须使用conda）:**

**推荐安装方式（所有系统）:**
```bash
conda install openbabel -c conda-forge
```

**⚠️ 错误方式（会导致功能缺失）:**
- ❌ `pip install openbabel` - 不存在此包
- ❌ `sudo apt-get install openbabel` - 缺少Python绑定
- ❌ `brew install openbabel` - 缺少Python绑定
- ❌ 从官网下载安装包 - 缺少Python绑定

**说明**: Open Babel的Python绑定只能通过conda安装，因为需要编译。

## 快速开始

### 导入方式

```python
# 标准导入（安装后使用）
import mol_conversion
from mol_conversion import FileConversion, MemoryConversion, GCNEncoding

# 快速便捷函数
from mol_conversion import quick_xyz_to_inchi, quick_split_multiframe
```

### 1. 文件转换模式（FileConversion）

适合需要将XYZ文件转换为其他格式文件的情况：

```python
from mol_conversion import FileConversion

# XYZ到SDF文件转换
FileConversion.xyz_to_sdf("molecule.xyz", "molecule.sdf")

# XYZ到PDB文件转换
FileConversion.xyz_to_pdb("molecule.xyz", "molecule.pdb")

# XYZ到MOL2文件转换
FileConversion.xyz_to_mol2("molecule.xyz", "molecule.mol2")

# XYZ到CIF文件转换
FileConversion.xyz_to_cif("molecule.xyz", "molecule.cif")

# XYZ文件到InChI字符串转换
inchi = FileConversion.xyz_to_inchi("molecule.xyz")
print(f"InChI: {inchi}")

# XYZ文件到SMILES字符串转换
smiles = FileConversion.xyz_to_smiles("molecule.xyz")
print(f"SMILES: {smiles}")

# 获取各种格式的字符串输出
pdb_str = FileConversion.xyz_to_pdb_string("molecule.xyz")
mol2_str = FileConversion.xyz_to_mol2_string("molecule.xyz")
cif_str = FileConversion.xyz_to_cif_string("molecule.xyz")

# 合并多个SDF文件
FileConversion.merge_sdf_files(
    sdf_list=["mol1.sdf", "mol2.sdf", "mol3.sdf"], 
    output_path="combined.sdf"
)
```

### 2. 内存转换模式（MemoryConversion）

适合需要链式转换和内存高效处理的场景：

```python
from memory_conversion import MemoryConversion

xyz_string = """14
Energy: -158.36316482
 C     1.965628    -0.044891     0.041336
 C    -1.928636     0.206074    -0.054869
 C    -0.547665    -0.397317     0.174957
 C     0.584659     0.558517    -0.188474
 H     0.479415     1.484663     0.402725
 H     0.479687     0.857511    -1.245797
 H    -0.442720    -0.696307     1.232280
 H    -0.442385    -1.323467    -0.416238
 H     2.101678    -0.324162     1.097821
 H     2.101862    -0.955927    -0.562083
 H     2.766586     0.658980    -0.226474
 H    -2.064698     0.485312    -1.111359
 H    -2.064872     1.117129     0.548524
 H    -2.729587    -0.497794     0.212971
"""

# XYZ字符串到MOL格式转换
mol_str = MemoryConversion.xyz_to_mol(xyz_string)
print(f"MOL格式: {mol_str[:100]}...")

# XYZ字符串到RDKit分子对象转换
rdkit_mol = MemoryConversion.xyz_to_rdkit(xyz_string)
print(f"RDKit分子对象: {type(rdkit_mol)}")

# XYZ字符串到键级矩阵转换
bond_matrix = MemoryConversion.xyz_to_bond_matrix(xyz_string)
print(f"键级矩阵形状: {bond_matrix.shape}")

# 完整的格式转换 - 文件格式字符串
file_formats = {
    'mol': MemoryConversion.xyz_to_mol_string(xyz_string),
    'sdf': MemoryConversion.xyz_to_sdf_string(xyz_string),
    'pdb': MemoryConversion.xyz_to_pdb_string(xyz_string),
    'mol2': MemoryConversion.xyz_to_mol2_string(xyz_string),
    'cif': MemoryConversion.xyz_to_cif_string(xyz_string)
}

for format_name, result in file_formats.items():
    print(f"{format_name.upper()}字符串长度: {len(result)}")

# 分子标识符转换
identifiers = {
    'inchi': MemoryConversion.xyz_to_inchi(xyz_string),
    'inchikey': MemoryConversion.xyz_to_inchikey(xyz_string),
    'smiles': MemoryConversion.xyz_to_smiles(xyz_string)
}

for name, value in identifiers.items():
    print(f"{name.upper()}: {value}")
```

### 3. 分子分析工具（GCN编码和标准化）

```python
from utils import GCNEncoding, XYZStandardizer

# 获取GCN编码
encoder = GCNEncoding(xyz_string)
encodings = encoder.encodings

print("GCN0编码:", encodings['gcn0'])
print("GCN1编码:", encodings['gcn1'])
print("GCN2编码:", encodings['gcn2'])

# 标准化分子坐标
standardizer = XYZStandardizer(xyz_string)
standardized_xyz = standardizer.to_standard_xyz()
print("标准化后的XYZ:")
print(standardized_xyz)
```

## 支持的转换格式

| 输入格式 | 输出格式 | FileConversion类 | MemoryConversion类 | 描述 |
|---------|---------|------------------|-------------------|------|
| XYZ文件 | SDF文件 | `xyz_to_sdf()` | - | 结构数据文件 |
| XYZ文件 | MOL文件 | `xyz_to_mol2()` | - | MOL2文件格式 |
| XYZ文件 | PDB文件 | `xyz_to_pdb()` | - | 蛋白质数据库格式 |
| XYZ文件 | CIF文件 | `xyz_to_cif()` | - | 晶体信息框架格式 |
| XYZ文件 | InChI | `xyz_to_inchi()` | - | IUPAC标准化学标识符 |
| XYZ文件 | InChIKey | `xyz_to_inchikey()` | - | InChI的哈希键 |
| XYZ文件 | SMILES | `xyz_to_smiles()` | - | 简化分子线性输入规范 |
| XYZ字符串 | MOL字符串 | - | `xyz_to_mol_string()` | MOL文件格式 |
| XYZ字符串 | SDF字符串 | - | `xyz_to_sdf_string()` | 结构数据文件 |
| XYZ字符串 | PDB字符串 | - | `xyz_to_pdb_string()` | 蛋白质数据库格式 |
| XYZ字符串 | MOL2字符串 | - | `xyz_to_mol2_string()` | MOL2文件格式 |
| XYZ字符串 | CIF字符串 | - | `xyz_to_cif_string()` | 晶体信息框架格式 |
| XYZ字符串 | InChI | - | `xyz_to_inchi()` | IUPAC标准化学标识符 |
| XYZ字符串 | InChIKey | - | `xyz_to_inchikey()` | InChI的哈希键 |
| XYZ字符串 | SMILES | - | `xyz_to_smiles()` | 简化分子线性输入规范 |
| XYZ字符串 | RDKit对象 | - | `xyz_to_rdkit()` | RDKit分子对象 |
| XYZ字符串 | 键级矩阵 | - | `xyz_to_bond_matrix()` | NumPy矩阵表示 |

## 选择指南

### 使用文件转换（mol_conversion.py）的情况：
- 需要输出文件到磁盘
- 处理大型分子文件
- 简单的单次转换
- 需要完整的Open Babel格式支持

### 使用内存转换（memory_conversion.py）的情况：
- 需要链式转换处理
- 内存充足且追求性能
- 批量处理中小型分子
- 需要访问中间结果（如键级矩阵）

### 使用分子分析工具（utils.py）的情况：
- 计算分子特征和编码
- 标准化分子坐标
- 图神经网络相关处理

## 性能考虑

- **内存转换**: 速度更快但内存占用较高
- **文件转换**: 内存占用低但涉及I/O开销
- **GCN编码**: 基于原子邻居关系，适合分子机器学习

## 示例文件

查看 `mol_conversion.ipynb` 获取完整的使用示例和可视化结果。

## 许可证

本项目遵循MIT许可证。

## 贡献

欢迎提交Issue和Pull Request来改进这个工具！