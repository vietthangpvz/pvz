# 1. Cài đặt thư viện (Chỉ cần chạy 1 lần duy nhất)
!pip install rdkit pandas

from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools, DataStructs
import pandas as pd

# 2. Giả lập dữ liệu đầu vào (SMILES)
# Giả sử ta có 4 chất, trong đó Aspirin bị lặp lại với 2 cách viết SMILES khác nhau
data = {
    'Compound_Name': ['Aspirin_1', 'Caffeine', 'Aspirin_2', 'Paracetamol'],
    'SMILES': [
        'CC(=O)OC1=CC=CC=C1C(=O)O',  # Aspirin cách 1
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', # Caffeine
        'O=C(C)Oc1ccccc1C(=O)O',      # Aspirin cách 2 (khác cách viết nhưng cùng 1 chất)
        'CC(=O)NC1=CC=C(O)C=C1'       # Paracetamol
    ]
}
df = pd.DataFrame(data)
print("--- Dữ liệu thô ban đầu ---")
print(df)

# 3. Sử dụng InChIKey để lọc trùng (Deduplication)
# Chuyển SMILES -> Mol -> InChIKey
df['Mol'] = df['SMILES'].apply(Chem.MolFromSmiles)
df['InChIKey'] = df['Mol'].apply(Chem.MolToInchiKey)

# Loại bỏ các dòng có InChIKey trùng nhau
df_clean = df.drop_duplicates(subset=['InChIKey']).copy()
print("\n--- Dữ liệu sau khi lọc trùng bằng InChIKey ---")
print(df_clean[['Compound_Name', 'InChIKey']])

# 4. Sử dụng Fingerprints để tìm chất giống Aspirin nhất
# Tạo Morgan Fingerprint cho tất cả các chất còn lại
df_clean['Fingerprint'] = df_clean['Mol'].apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=1024))

# Lấy Fingerprint của Aspirin làm chuẩn (mốc so sánh)
aspirin_fp = df_clean.iloc[0]['Fingerprint']

# Tính độ tương đồng Tanimoto của các chất khác so với Aspirin
df_clean['Similarity_to_Aspirin'] = df_clean['Fingerprint'].apply(lambda x: DataStructs.TanimotoSimilarity(aspirin_fp, x))

# 5. Hiển thị kết quả đẹp mắt với PandasTools
print("\n--- Kết quả phân tích cuối cùng ---")
# Thêm cột hình ảnh để xem trực quan trên bảng
PandasTools.AddMoleculeColumnToFrame(df_clean, 'SMILES', 'Structure')
df_clean[['Compound_Name', 'Similarity_to_Aspirin', 'Structure']]
