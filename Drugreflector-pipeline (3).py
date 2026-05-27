import scanpy as sc
import pandas as pd
import numpy as np
import os
import json
import drugreflector as dr 
def prepare_drugreflector_input(
    adata_path,
    group_col,
    group1,       # 正常/对照组
    group2,       # 疾病/实验组
    output_dir = "./data/drugreflector_input/",
    cache_file = "gene_rename_cache.json"
):
    """
    从单细胞数据准备 DrugReflector 标准输入文件
    
    Parameters
    ----------
    adata_path : str
        h5ad 文件路径
    group_col : str
        分组列名（如 'cell_type' 或 'condition'）
    group1 : str
        对照组标签
    group2 : str
        实验组标签（疾病状态）
    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("===== 准备 DrugReflector 输入文件 =====")
    
    # 1. 加载数据
    print(f"\n[1/5] 加载数据: {adata_path}")
    adata = sc.read(adata_path)
    print(f"      数据形状: {adata.shape}")
    
    # 2. 自动更新基因名
    print("\n[2/5] 更新基因名至 HGNC 标准")
    
    def auto_rename_genes(var_names, cache_file):
        names = list(var_names)
        if os.path.exists(cache_file):
            with open(cache_file) as f:
                rename = json.load(f)
            print(f"      从缓存读取 {len(rename)} 条映射")
        else:
            import mygene
            mg = mygene.MyGeneInfo()
            result = mg.querymany(names, scopes='symbol,alias',
                                  fields='symbol', species='human',
                                  verbose=False)
            rename = {r['query']: r['symbol'] 
                      for r in result
                      if 'symbol' in r and r['symbol'] != r['query']}
            with open(cache_file, 'w') as f:
                json.dump(rename, f, indent=2)
            print(f"      联网查询完成，更新 {len(rename)} 个基因名")
        return [rename.get(g, g) for g in names]
    
    adata.var_names = auto_rename_genes(adata.var_names, cache_file)
    
    # 3. 标准化
    print("\n[3/5] 数据标准化")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print("      normalize_total(1e4) + log1p 完成")
    
    # 4. 计算 v-score
    print(f"\n[4/5] 计算 V-score: {group1} → {group2}")
    
    mask1 = adata.obs[group_col] == group1
    mask2 = adata.obs[group_col] == group2
    
    if mask1.sum() == 0:
        raise ValueError(f"找不到 {group1}，请检查 group_col='{group_col}'")
    if mask2.sum() == 0:
        raise ValueError(f"找不到 {group2}，请检查 group_col='{group_col}'")
    
    print(f"      {group1}: {mask1.sum()} 细胞")
    print(f"      {group2}: {mask2.sum()} 细胞")
    
    expr1 = adata[mask1].X
    expr2 = adata[mask2].X
    if hasattr(expr1, 'toarray'):
        expr1 = expr1.toarray()
        expr2 = expr2.toarray()
    
    vscores = pd.Series(
        expr2.mean(axis=0) - expr1.mean(axis=0),
        index = adata.var_names,
        name  = f"{group1}->{group2}"
    )
    
    print(f"      V-score 计算完成，{len(vscores)} 个基因")
    print(f"      Top 5 上调: {list(vscores.nlargest(5).index)}")
    print(f"      Top 5 下调: {list(vscores.nsmallest(5).index)}")
    
    # 5. 保存文件
    print("\n[5/5] 保存输入文件")
    
    vscore_path = os.path.join(output_dir, "vscores.csv")
    vscores.to_csv(vscore_path, header=True)
    print(f"      V-score 已保存: {vscore_path}")
    
    # 同时保存 AnnData 格式
    vscore_adata = sc.AnnData(
        X   = vscores.values.reshape(1, -1),
        var = pd.DataFrame(index=vscores.index),
        obs = pd.DataFrame({'comparison': [vscores.name]})
    )
    vscore_adata.write(os.path.join(output_dir, "vscores.h5ad"))
    
    print(f"\n===== 输入文件准备完成 =====")
    print(f"输出目录: {output_dir}")
    
    return vscores

vscores = prepare_drugreflector_input(
    adata_path = "./drugreflector/my_expression_data.h5ad",
    group_col  = "cellType",     # 根据实际数据调整
    group1     = "Breast_Cancer_Cells",       # 原发灶乳腺癌细胞
    group2     = "T_cells",     # T
)

def check_data_quality(vscores, model):
    """运行前数据质量检查"""
    
    print("===== 数据质量检查 =====")
    
    # 1. 基因覆盖率
    coverage = model.check_gene_coverage(vscores.index)
    print(f"\n基因覆盖率: {coverage['coverage_percent']:.1f}%")
    print(f"  模型需要: 978 个 landmark 基因")
    print(f"  数据包含: {coverage['total_found']} 个")
    print(f"  缺失:     {coverage['total_input'] - coverage['total_found']} 个")
    
    if coverage['coverage_percent'] < 80:
        print("\n⚠️  警告：基因覆盖率低于80%，预测结果可信度降低")
        print("   建议检查基因命名格式是否为 HGNC 标准符号")
    else:
        print("\n✓ 基因覆盖率良好，预测结果可信")
    
    # 2. V-score 分布检查
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    vscores.hist(bins=50, ax=axes[0], color='#4DBBD5', edgecolor='white')
    axes[0].set_title('V-score 分布')
    axes[0].set_xlabel('V-score')
    axes[0].axvline(0, color='red', linestyle='--')
    
    import scipy.stats as stats
    stats.probplot(vscores, plot=axes[1])
    axes[1].set_title('Q-Q 图（正态性检验）')
    
    plt.tight_layout()
    plt.savefig('./figures/vscore_distribution.pdf', bbox_inches='tight')
    
    print(f"\nV-score 统计:")
    print(f"  均值:   {vscores.mean():.4f}")
    print(f"  标准差: {vscores.std():.4f}")
    print(f"  范围:   [{vscores.min():.4f}, {vscores.max():.4f}]")
    
    return coverage

def check_data_quality(vscores, model):
    """运行前数据质量检查"""
    
    print("===== 数据质量检查 =====")
    
    # 1. 基因覆盖率
    coverage = model.check_gene_coverage(vscores.index)
    print(f"\n基因覆盖率: {coverage['coverage_percent']:.1f}%")
    print(f"  模型需要: 978 个 landmark 基因")
    print(f"  数据包含: {coverage['total_found']} 个")
    print(f"  缺失:     {coverage['total_input'] - coverage['total_found']} 个")
    
    if coverage['coverage_percent'] < 80:
        print("\n  警告：基因覆盖率低于80%，预测结果可信度降低")
        print("   建议检查基因命名格式是否为 HGNC 标准符号")
    else:
        print("\n✓ 基因覆盖率良好，预测结果可信")
    
    # 2. V-score 分布检查
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    vscores.hist(bins=50, ax=axes[0], color='#4DBBD5', edgecolor='white')
    axes[0].set_title('V-score 分布')
    axes[0].set_xlabel('V-score')
    axes[0].axvline(0, color='red', linestyle='--')
    
    import scipy.stats as stats
    stats.probplot(vscores, plot=axes[1])
    axes[1].set_title('Q-Q 图（正态性检验）')
    
    plt.tight_layout()
    plt.savefig('./figures/vscore_distribution.pdf', bbox_inches='tight')
    
    print(f"\nV-score 统计:")
    print(f"  均值:   {vscores.mean():.4f}")
    print(f"  标准差: {vscores.std():.4f}")
    print(f"  范围:   [{vscores.min():.4f}, {vscores.max():.4f}]")
    
    return coverage

import sys, types, torch, warnings, json, os
import torch.nn as nn

warnings.filterwarnings('ignore', message='.*InconsistentVersionWarning.*')
warnings.filterwarnings('ignore', category=UserWarning, module='sklearn')
warnings.filterwarnings('ignore', message='.*missing.*genes.*')

class _CallableModule(types.ModuleType):
    def __call__(self, *args, **kwargs):
        return None
    def __getattr__(self, name):
        full_name = self.__name__ + '.' + name
        if full_name not in sys.modules:
            sub = _CallableModule(full_name)
            sys.modules[full_name] = sub
            object.__setattr__(self, name, sub)
        return sys.modules[full_name]

class _MLP(nn.Module):
    def __init__(self, *args, **kwargs):
        super().__init__()
        self.net = nn.Sequential()
    def forward(self, x): return self.net(x)

class _FocalLoss(nn.Module):
    def __init__(self, *args, **kwargs): super().__init__()
    def forward(self, x, y=None): return x

cifra = _CallableModule('cifra')
sys.modules['cifra'] = cifra
for _mod in ['cifra.torch','cifra.torch.losses',
             'cifra.torch.losses.focal_loss',
             'cifra.torch.models','cifra.torch.utils',
             'cifra.utils','cifra.data']:
    _m = _CallableModule(_mod)
    _m.FocalLoss = _FocalLoss
    _m.MLP = _MLP
    sys.modules[_mod] = _m

# ===== 正式代码 =====
import numpy as np
import pandas as pd
import scanpy as sc
import drugreflector as dr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def run_drugreflector_screening(
    vscores,
    checkpoint_dir  = "./checkpoints/",
    output_dir      = "./results/drugreflector/",
    n_top           = 200,
    study_name      = "Breast_Cancer_Brain_Met"
):
    """
    完整 DrugReflector 筛选流程
    
    Parameters
    ----------
    vscores : pd.Series
        V-score 向量，index 为基因名
    checkpoint_dir : str
        模型权重目录
    output_dir : str
        结果输出目录
    n_top : int
        输出 Top N 化合物
    study_name : str
        研究名称（用于文件命名）
    """
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs("./figures/", exist_ok=True)
    
    print(f"===== DrugReflector 药物筛选: {study_name} =====")
    
    # 1. 加载模型
    print("\n[1/4] 加载模型权重...")
    drugreflector/checkpoints
    checkpoint_dir = "./drugreflector/checkpoints"
    model_paths = [
    os.path.join(checkpoint_dir, 'model_fold_0.pt'),
    os.path.join(checkpoint_dir, 'model_fold_1.pt'),
    os.path.join(checkpoint_dir, "model_fold_2.pt")
    # ... 其他模型路径
]
    model = dr.DrugReflector(checkpoint_paths=model_paths)
    print(f"      模型加载成功")
    print(f"      输入基因: {len(model.model.dimensions['var_names'][0])}")
    print(f"      输出化合物: {len(model.model.dimensions['output_names'])}")
    
    # 2. 数据质量检查
    print("\n[2/4] 检查数据质量...")
    coverage = model.check_gene_coverage(vscores.index)
    print(f"      基因覆盖率: {coverage['coverage_percent']:.1f}%")
    
    # 3. 预测
    print(f"\n[3/4] 运行预测（Top {n_top} 化合物）...")
    predictions = model.predict(vscores, n_top=n_top)
    
    # 整理结果
    rank_col  = [c for c in predictions.columns if c[0] == 'rank'][0]
    prob_col  = [c for c in predictions.columns if c[0] == 'prob'][0]
    logit_col = [c for c in predictions.columns if c[0] == 'logit'][0]
    
    results_df = pd.DataFrame({
        'compound_id': predictions.index,
        'rank':        predictions[rank_col].values,
        'probability': predictions[prob_col].values,
        'logit_score': predictions[logit_col].values,
    }).sort_values('rank').reset_index(drop=True)
    
    results_df['rank_percentile'] = (
        1 - results_df['rank'] / len(model.model.dimensions['output_names'])
    ) * 100
    
    print(f"\n      Top 10 候选化合物:")
    print(f"      {'排名':<6} {'化合物ID':<20} {'概率':<10} {'百分位':<10}")
    print(f"      {'-'*50}")
    for _, row in results_df.head(10).iterrows():
        print(f"      {int(row['rank']):<6} {row['compound_id']:<20} "
              f"{row['probability']:.4f}    "
              f"Top {100-row['rank_percentile']:.1f}%")
    
    # 4. 保存结果
    print(f"\n[4/4] 保存结果...")
    
    out_csv = os.path.join(output_dir, f"{study_name}_predictions.csv")
    results_df.to_csv(out_csv, index=False)
    print(f"      CSV: {out_csv}")
    
    # 保存完整预测（所有9597化合物）
    full_predictions = model.predict(vscores, n_top=None)
    full_rank_col = [c for c in full_predictions.columns if c[0]=='rank'][0]
    full_results = pd.DataFrame({
        'compound_id': full_predictions.index,
        'rank':        full_predictions[full_rank_col].values,
        'probability': full_predictions[[c for c in full_predictions.columns 
                                         if c[0]=='prob'][0]].values,
    }).sort_values('rank').reset_index(drop=True)
    
    full_out = os.path.join(output_dir, f"{study_name}_all_compounds.csv")
    full_results.to_csv(full_out, index=False)
    print(f"      全量: {full_out}")
    
    return results_df, full_results


# 运行筛选
vscores = pd.read_csv(
    "./data/drugreflector_input/vscores.csv", index_col=0
).squeeze()

top_results, all_results = run_drugreflector_screening(
    vscores     = vscores,
    n_top       = 200,
    study_name  = "BC"   # Breast Cancer Brain Metastasis
)


# ===== 正式代码 =====
import numpy as np
import pandas as pd
import scanpy as sc
import drugreflector as dr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import json

def run_drugreflector_screening(
    vscores,
    checkpoint_dir  = "./checkpoints/",
    output_dir      = "./results/drugreflector/",
    n_top           = 200,
    study_name      = "Breast_Cancer_Brain_Met"
):
    """
    完整 DrugReflector 筛选流程
    
    Parameters
    ----------
    vscores : pd.Series
        V-score 向量，index 为基因名
    checkpoint_dir : str
        模型权重目录
    output_dir : str
        结果输出目录
    n_top : int
        输出 Top N 化合物
    study_name : str
        研究名称（用于文件命名）
    """
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs("./figures/", exist_ok=True)
    
    print(f"===== DrugReflector 药物筛选: {study_name} =====")
    
    # 1. 加载模型
    print("\n[1/4] 加载模型权重...")
    model = dr.DrugReflector(checkpoint_paths=[
        os.path.join(checkpoint_dir, "model_fold_0.pt"),
        os.path.join(checkpoint_dir, "model_fold_1.pt"),
        os.path.join(checkpoint_dir, "model_fold_2.pt")
    ])
    print(f"      模型加载成功")
    print(f"      输入基因: {len(model.model.dimensions['var_names'][0])}")
    print(f"      输出化合物: {len(model.model.dimensions['output_names'])}")
    
    # 2. 数据质量检查
    print("\n[2/4] 检查数据质量...")
    coverage = model.check_gene_coverage(vscores.index)
    print(f"      基因覆盖率: {coverage['coverage_percent']:.1f}%")
    
    # 3. 预测
    print(f"\n[3/4] 运行预测（Top {n_top} 化合物）...")
    predictions = model.predict(vscores, n_top=n_top)
    
    # 整理结果
    rank_col  = [c for c in predictions.columns if c[0] == 'rank'][0]
    prob_col  = [c for c in predictions.columns if c[0] == 'prob'][0]
    logit_col = [c for c in predictions.columns if c[0] == 'logit'][0]
    
    results_df = pd.DataFrame({
        'compound_id': predictions.index,
        'rank':        predictions[rank_col].values,
        'probability': predictions[prob_col].values,
        'logit_score': predictions[logit_col].values,
    }).sort_values('rank').reset_index(drop=True)
    
    results_df['rank_percentile'] = (
        1 - results_df['rank'] / len(model.model.dimensions['output_names'])
    ) * 100
    
    print(f"\n      Top 10 候选化合物:")
    print(f"      {'排名':<6} {'化合物ID':<20} {'概率':<10} {'百分位':<10}")
    print(f"      {'-'*50}")
    for _, row in results_df.head(10).iterrows():
        print(f"      {int(row['rank']):<6} {row['compound_id']:<20} "
              f"{row['probability']:.4f}    "
              f"Top {100-row['rank_percentile']:.1f}%")
    
    # 4. 保存结果
    print(f"\n[4/4] 保存结果...")
    
    out_csv = os.path.join(output_dir, f"{study_name}_predictions.csv")
    results_df.to_csv(out_csv, index=False)
    print(f"      CSV: {out_csv}")
    
    # 保存完整预测（所有9597化合物）
    full_predictions = model.predict(vscores, n_top=None)
    full_rank_col = [c for c in full_predictions.columns if c[0]=='rank'][0]
    full_results = pd.DataFrame({
        'compound_id': full_predictions.index,
        'rank':        full_predictions[full_rank_col].values,
        'probability': full_predictions[[c for c in full_predictions.columns 
                                         if c[0]=='prob'][0]].values,
    }).sort_values('rank').reset_index(drop=True)
    
    full_out = os.path.join(output_dir, f"{study_name}_all_compounds.csv")
    full_results.to_csv(full_out, index=False)
    print(f"      全量: {full_out}")
    
    return results_df, full_results


# 运行筛选
vscores = pd.read_csv(
    "./data/drugreflector_input/vscores.csv", index_col=0
).squeeze()

top_results, all_results = run_drugreflector_screening(
    vscores     = vscores,
    n_top       = 200,
    study_name  = "BC"   # Breast Cancer Brain Metastasis
)



import scanpy as sc
import pandas as pd
import numpy as np
import os
import json
import matplotlib.pyplot as plt
import scipy.stats as stats 
# ================================================================
# DrugReflector 完整运行脚本（生产版本）
# ================================================================

# ===== cifra 模块动态注入 =====
import sys, types, torch, warnings, json, os
import torch.nn as nn

warnings.filterwarnings('ignore', message='.*InconsistentVersionWarning.*')
warnings.filterwarnings('ignore', category=UserWarning, module='sklearn')
warnings.filterwarnings('ignore', message='.*missing.*genes.*')

class _CallableModule(types.ModuleType):
    def __call__(self, *args, **kwargs):
        return None
    def __getattr__(self, name):
        full_name = self.__name__ + '.' + name
        if full_name not in sys.modules:
            sub = _CallableModule(full_name)
            sys.modules[full_name] = sub
            object.__setattr__(self, name, sub)
        return sys.modules[full_name]

class _MLP(nn.Module):
    def __init__(self, *args, **kwargs):
        super().__init__()
        self.net = nn.Sequential()
    def forward(self, x): return self.net(x)

class _FocalLoss(nn.Module):
    def __init__(self, *args, **kwargs): super().__init__()
    def forward(self, x, y=None): return x

cifra = _CallableModule('cifra')
sys.modules['cifra'] = cifra
for _mod in ['cifra.torch','cifra.torch.losses',
             'cifra.torch.losses.focal_loss',
             'cifra.torch.models','cifra.torch.utils',
             'cifra.utils','cifra.data']:
    _m = _CallableModule(_mod)
    _m.FocalLoss = _FocalLoss
    _m.MLP = _MLP
    sys.modules[_mod] = _m
# ==============================================================================
# 1. 函数定义部分
# ==============================================================================

def prepare_drugreflector_input(
    adata_path,
    group_col,
    group1,       # 正常/对照组
    group2,       # 疾病/实验组
    output_dir = "./data/drugreflector_input/",
    cache_file = "gene_rename_cache.json"
):
    """
    从单细胞数据准备 DrugReflector 标准输入文件
    
    Parameters
    ----------
    adata_path : str
        h5ad 文件路径
    group_col : str
        分组列名（如 'cell_type' 或 'condition'）
    group1 : str
        对照组标签
    group2 : str
        实验组标签（疾病状态）
    output_dir : str
        输出目录
    cache_file : str
        基因名转换的缓存文件路径
    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("===== 准备 DrugReflector 输入文件 =====")
    
    # 1. 加载数据
    print(f"\n[1/5] 加载数据: {adata_path}")
    adata = sc.read(adata_path)
    print(f"      数据形状: {adata.shape}")
    
    # 2. 自动更新基因名
    print("\n[2/5] 更新基因名至 HGNC 标准")
    
    def auto_rename_genes(var_names, cache_file):
        names = list(var_names)
        if os.path.exists(cache_file):
            with open(cache_file) as f:
                rename = json.load(f)
            print(f"      从缓存读取 {len(rename)} 条映射")
        else:
            try:
                import mygene
                mg = mygene.MyGeneInfo()
                result = mg.querymany(names, scopes='symbol,alias',
                                      fields='symbol', species='human',
                                      verbose=False)
                rename = {r['query']: r['symbol'] 
                          for r in result
                          if 'symbol' in r and r['symbol'] != r['query']}
                with open(cache_file, 'w') as f:
                    json.dump(rename, f, indent=2)
                print(f"      联网查询完成，更新 {len(rename)} 个基因名")
            except ImportError:
                print("      'mygene' 未安装，跳过基因名更新。请运行 `pip install mygene`")
                rename = {}
        
        # 使用 pd.Index.map 更高效地进行重命名和去重
        new_var_names = pd.Index(names).map(lambda g: rename.get(g, g))
        adata.var_names_make_unique() # 在重命名后处理重复基因名
        return new_var_names

    adata.var_names = auto_rename_genes(adata.var_names, cache_file)
    adata.var_names_make_unique() # 确保最终没有重复的基因名
    
    # 3. 标准化
    print("\n[3/5] 数据标准化")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print("      normalize_total(1e4) + log1p 完成")
    
    # 4. 计算 v-score
    print(f"\n[4/5] 计算 V-score: {group1} vs {group2}")
    
    mask1 = adata.obs[group_col] == group1
    mask2 = adata.obs[group_col] == group2
    
    if mask1.sum() == 0:
        raise ValueError(f"在列 '{group_col}' 中找不到分组 '{group1}'，请检查分组名称是否正确。")
    if mask2.sum() == 0:
        raise ValueError(f"在列 '{group_col}' 中找不到分组 '{group2}'，请检查分组名称是否正确。")
    
    print(f"      {group1} (对照组): {mask1.sum()} 细胞")
    print(f"      {group2} (实验组): {mask2.sum()} 细胞")
    
    # 使用 scanpy 的 rank_genes_groups 计算更稳健，但为了保持与原版一致，这里仍用均值差
    # 注意：v-score 定义是 实验组 - 对照组
    expr1 = adata[mask1].X
    expr2 = adata[mask2].X
    if hasattr(expr1, 'toarray'):
        expr1 = expr1.toarray()
        expr2 = expr2.toarray()
    
    vscores = pd.Series(
        expr2.mean(axis=0) - expr1.mean(axis=0),
        index = adata.var_names,
        name  = f"{group2}_vs_{group1}"
    )
    
    print(f"      V-score 计算完成，共 {len(vscores)} 个基因")
    print(f"      Top 5 上调 (在 {group2} 中高表达): {list(vscores.nlargest(5).index)}")
    print(f"      Top 5 下调 (在 {group1} 中高表达): {list(vscores.nsmallest(5).index)}")
    
    # 5. 保存文件
    print("\n[5/5] 保存输入文件")
    
    vscore_path = os.path.join(output_dir, "vscores.csv")
    vscores.to_csv(vscore_path, header=True)
    print(f"      V-score 已保存: {vscore_path}")
    
    print(f"\n===== 输入文件准备完成 =====")
    print(f"输出目录: {output_dir}")
    
    return vscores

def check_data_quality(vscores, model, output_dir="./figures"):
    """
    对生成的 v-score 进行质量检查，并生成可视化报告。

    Parameters
    ----------
    vscores : pd.Series
        计算得到的 v-score。
    model : drugreflector.Model
        已加载的 DrugReflector 模型对象。
    output_dir : str
        保存可视化图表的目录。
    """
    
    print("\n===== 数据质量检查 =====")
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 基因覆盖率
    # 注意：新版 drugreflector 可能使用 model.gsea.check_gene_coverage
    # 这里我们假设是 model.check_gene_coverage
    try:
        coverage = model.check_gene_coverage(vscores.index)
        print(f"\n[QC 1/2] 基因覆盖率检查")
        print(f"  模型需要: {coverage['total_input']} 个 landmark 基因")
        print(f"  您的数据提供了: {coverage['total_found']} 个")
        print(f"  覆盖率: {coverage['coverage_percent']:.1f}%")
        
        if coverage['coverage_percent'] < 80:
            print("\n  ⚠️  警告：基因覆盖率低于80%，预测结果可信度会降低。")
            print("     请检查基因名是否为最新的 HGNC 标准符号。")
        else:
            print("\n  ✓ 基因覆盖率良好，预测结果可信度较高。")
    except Exception as e:
        print(f"\n[QC 1/2] 基因覆盖率检查失败: {e}")
        print("  跳过此项检查。")


    # 2. V-score 分布检查与可视化
    print("\n[QC 2/2] V-score 分布检查与可视化")
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    # 直方图
    vscores.hist(bins=100, ax=axes[0], color='#4DBBD5', edgecolor='white', alpha=0.8)
    axes[0].set_title('V-score Distribution')
    axes[0].set_xlabel('V-score (log-fold change approximation)')
    axes[0].set_ylabel('Frequency')
    axes[0].axvline(0, color='red', linestyle='--')
    
    # Q-Q图
    stats.probplot(vscores, dist="norm", plot=axes[1])
    axes[1].get_lines()[0].set_markerfacecolor('#4DBBD5')
    axes[1].get_lines()[0].set_markeredgecolor('#4DBBD5')
    axes[1].get_lines()[1].set_color('red')
    axes[1].set_title('Normal Q-Q Plot')
    
    plt.tight_layout()
    
    # 保存 PDF 文件
    pdf_path = os.path.join(output_dir, 'vscore_distribution_qc.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"  可视化报告已保存: {pdf_path}")
    
    # 终端打印统计信息
    print(f"\n  V-score 统计摘要:")
    print(f"    均值 (Mean):   {vscores.mean():.4f}")
    print(f"    标准差 (Std): {vscores.std():.4f}")
    print(f"    中位数 (Median): {vscores.median():.4f}")
    print(f"    范围 (Range):   [{vscores.min():.4f}, {vscores.max():.4f}]")
    
    print("\n===== 质量检查完成 =====")

# CORRECTED DRUG SCREENING FUNCTION
# Please use this version in your main script.
# ==============================================================================
def run_drugreflector_screening(
    vscores,
    analysis_name="screening_results",
    checkpoint_dir="./drugreflector/checkpoints",
    output_dir="./data/drugreflector_output/"
):
    """
    使用 v-score 运行 DrugReflector 药物筛选。
    """
    print(f"\n===== DrugReflector 药物筛选: {analysis_name} =====")
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 加载模型权重
    print("\n[1/4] 加载模型权重...")
    model_paths = [
        os.path.join(checkpoint_dir, 'model_fold_0.pt'),
        os.path.join(checkpoint_dir, 'model_fold_1.pt'),
        os.path.join(checkpoint_dir, 'model_fold_2.pt')
    ]
    model = dr.DrugReflector(checkpoint_paths=model_paths)
    print("      模型加载成功!")

    # 2. Run data quality check
    qc_output_dir = os.path.join(output_dir, f"{analysis_name}_qc_reports")
    check_data_quality(vscores, model, output_dir=qc_output_dir)

    # 3. Perform drug prediction
    print("\n[3/4] 进行药物预测...")
    results_df = model.predict(vscores)
    print("      预测完成！")

    # 4. Save results
    print("\n[4/4] 保存预测结果...")
    
    # --- The definitive fix for MultiIndex columns ---
    
    # The output has MultiIndex columns. We want to sort by 'logit'.
    # The sample name is the second level of the index.
    sample_name = results_df.columns.levels[1][0]
    score_column_tuple = ('logit', sample_name)
    
    print(f"      使用列 '{score_column_tuple}' 作为分数进行排序。")

    # Sort the original DataFrame by the logit score (lower is better)
    sorted_df = results_df.sort_values(by=score_column_tuple, ascending=True)

    # For simplicity in the output CSV, let's flatten the MultiIndex columns.
    # We will join the two levels of the column names with an underscore.
    # e.g., ('logit', 'sample1') becomes 'logit_sample1'
    final_results = sorted_df.copy()
    final_results.columns = ['_'.join(col).strip() for col in final_results.columns.values]
    
    # Let's rename the primary score column to just 'score' for clarity
    final_results.rename(columns={f'logit_{sample_name}': 'score'}, inplace=True)
    
    # --- End of the fix ---

    all_results_path = os.path.join(output_dir, f"{analysis_name}_all_results.csv")
    final_results.to_csv(all_results_path, index=True)
    print(f"      所有药物的预测结果已保存: {all_results_path}")
    
    top30_results_path = os.path.join(output_dir, f"{analysis_name}_top30_results.csv")
    final_results.head(30).to_csv(top30_results_path, index=True)
    print(f"      Top 30 药物预测结果已保存: {top30_results_path}")
    
    print("\n===== 筛选完成 =====")
    
    return final_results.head(30), final_results

# ==============================================================================
# AND ENSURE YOUR MAIN CALL IS CORRECT
# ==============================================================================
# vscores = prepare_drugreflector_input(...) # Make sure this runs first

top_results, all_results = run_drugreflector_screening(
    vscores=vscores,
    analysis_name="BC", # Or whatever you named it
    checkpoint_dir="./drugreflector/checkpoints" # Pass the path as a string argument
)

print(top_results.head())
# ==============================================================================
# 2. 主程序执行部分
# ==============================================================================

if __name__ == '__main__':
    
    # --- 步骤 1: 准备 v-score ---
    # 调用第一个函数，生成 v-score
    vscores = prepare_drugreflector_input(
        adata_path = "./drugreflector/my_expression_data.h5ad", # 确保此路径正确
        group_col  = "cellType",     # 确保此列名在你的 .obs 中存在
        group1     = "Breast_Cancer_Cells",       # 对照组
        group2     = "T_cells",     # 实验组
    )

    # --- 步骤 2: 加载 DrugReflector 模型 ---
    # 这是调用质量检查前必须的一步
    
 # 1. 加载模型
   
    print("\n[1/4] 加载模型权重...")
    checkpoint_dir = "./drugreflector/checkpoints"
    model_paths = [
    os.path.join(checkpoint_dir, 'model_fold_0.pt'),
    os.path.join(checkpoint_dir, 'model_fold_1.pt'),
    os.path.join(checkpoint_dir, "model_fold_2.pt")
    # ... 其他模型路径
]
    model = dr.DrugReflector(checkpoint_paths=model_paths)
    print(f"      模型加载成功")
    print(f"      输入基因: {len(model.model.dimensions['var_names'][0])}")
    print(f"      输出化合物: {len(model.model.dimensions['output_names'])}")
    # --- 步骤 3: 执行质量检查 ---
    # 关键一步：调用 check_data_quality 函数，并传入上一步生成的 vscores 和 model
    check_data_quality(
        vscores = vscores, 
        model   = model,
        output_dir = "./figures" # 指定图表输出目录
    )
try:
    import drugreflector as dr
except ImportError:
    print("错误：'drugreflector' 库未安装或找不到。")
    print("请按照其 GitHub 页面的指引进行安装。")、

# ===== 正式代码 =====
import numpy as np
import pandas as pd
import scanpy as sc
import drugreflector as dr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def run_drugreflector_screening(
    vscores,
    analysis_name="screening_results",
    checkpoint_dir="./drugreflector/checkpoints",
    output_dir="./data/drugreflector_output/"
):
    """
    使用 v-score 运行 DrugReflector 药物筛选。
    """
    print(f"\n===== DrugReflector 药物筛选: {analysis_name} =====")
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 加载模型权重
    print("\n[1/4] 加载模型权重...")
    model_paths = [
        os.path.join(checkpoint_dir, 'model_fold_0.pt'),
        os.path.join(checkpoint_dir, 'model_fold_1.pt'),
        os.path.join(checkpoint_dir, 'model_fold_2.pt')
    ]
    model = dr.DrugReflector(checkpoint_paths=model_paths)
    print("      模型加载成功!")

    # 2. Run data quality check
    qc_output_dir = os.path.join(output_dir, f"{analysis_name}_qc_reports")
    check_data_quality(vscores, model, output_dir=qc_output_dir)

    # 3. Perform drug prediction
    print("\n[3/4] 进行药物预测...")
    results_df = model.predict(vscores)
    print("      预测完成！")

    # 4. Save results
    print("\n[4/4] 保存预测结果...")
    
    # --- The definitive fix for MultiIndex columns ---
    
    # The output has MultiIndex columns. We want to sort by 'logit'.
    # The sample name is the second level of the index.
    sample_name = results_df.columns.levels[1][0]
    score_column_tuple = ('logit', sample_name)
    
    print(f"      使用列 '{score_column_tuple}' 作为分数进行排序。")

    # Sort the original DataFrame by the logit score (lower is better)
    sorted_df = results_df.sort_values(by=score_column_tuple, ascending=True)

    # For simplicity in the output CSV, let's flatten the MultiIndex columns.
    # We will join the two levels of the column names with an underscore.
    # e.g., ('logit', 'sample1') becomes 'logit_sample1'
    final_results = sorted_df.copy()
    final_results.columns = ['_'.join(col).strip() for col in final_results.columns.values]
    
    # Let's rename the primary score column to just 'score' for clarity
    final_results.rename(columns={f'logit_{sample_name}': 'score'}, inplace=True)
    
    # --- End of the fix ---

    all_results_path = os.path.join(output_dir, f"{analysis_name}_all_results.csv")
    final_results.to_csv(all_results_path, index=True)
    print(f"      所有药物的预测结果已保存: {all_results_path}")
    
    top30_results_path = os.path.join(output_dir, f"{analysis_name}_top30_results.csv")
    final_results.head(30).to_csv(top30_results_path, index=True)
    print(f"      Top 30 药物预测结果已保存: {top30_results_path}")
    
    print("\n===== 筛选完成 =====")
    
    return final_results.head(30), final_results

def run_drugreflector_screening(
    vscores,
    analysis_name="screening_results",
    checkpoint_dir="./drugreflector/checkpoints",
    output_dir="./data/drugreflector_output/"
):
    """
    使用 v-score 运行 DrugReflector 药物筛选。
    """
    print(f"\n===== DrugReflector 药物筛选: {analysis_name} =====")
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 加载模型权重
    print("\n[1/4] 加载模型权重...")
    model_paths = [
        os.path.join(checkpoint_dir, 'model_fold_0.pt'),
        os.path.join(checkpoint_dir, 'model_fold_1.pt'),
        os.path.join(checkpoint_dir, 'model_fold_2.pt')
    ]
    model = dr.DrugReflector(checkpoint_paths=model_paths)
    print("      模型加载成功!")

    # 2. Run data quality check
    qc_output_dir = os.path.join(output_dir, f"{analysis_name}_qc_reports")
    check_data_quality(vscores, model, output_dir=qc_output_dir)

    # 3. Perform drug prediction
    print("\n[3/4] 进行药物预测...")
    results_df = model.predict(vscores)
    print("      预测完成！")

    # 4. Save results
    print("\n[4/4] 保存预测结果...")
    
    # --- The definitive fix for MultiIndex columns ---
    
    # The output has MultiIndex columns. We want to sort by 'logit'.
    # The sample name is the second level of the index.
    sample_name = results_df.columns.levels[1][0]
    score_column_tuple = ('logit', sample_name)
    
    print(f"      使用列 '{score_column_tuple}' 作为分数进行排序。")

    # Sort the original DataFrame by the logit score (lower is better)
    sorted_df = results_df.sort_values(by=score_column_tuple, ascending=True)

    # For simplicity in the output CSV, let's flatten the MultiIndex columns.
    # We will join the two levels of the column names with an underscore.
    # e.g., ('logit', 'sample1') becomes 'logit_sample1'
    final_results = sorted_df.copy()
    final_results.columns = ['_'.join(col).strip() for col in final_results.columns.values]
    
    # Let's rename the primary score column to just 'score' for clarity
    final_results.rename(columns={f'logit_{sample_name}': 'score'}, inplace=True)
    
    # --- End of the fix ---

    all_results_path = os.path.join(output_dir, f"{analysis_name}_all_results.csv")
    final_results.to_csv(all_results_path, index=True)
    print(f"      所有药物的预测结果已保存: {all_results_path}")
    
    top30_results_path = os.path.join(output_dir, f"{analysis_name}_top30_results.csv")
    final_results.head(30).to_csv(top30_results_path, index=True)
    print(f"      Top 30 药物预测结果已保存: {top30_results_path}")
    
    print("\n===== 筛选完成 =====")
    
    return final_results.head(30), final_results
# 运行筛选

vscores = pd.read_csv(
 "./data/drugreflector_input/vscores.csv", index_col=0
).squeeze()

# Call the function with the correct parameter names
top_results, all_results = run_drugreflector_screening(
    vscores = vscores,
    analysis_name = "BC"  # Use 'analysis_name' instead of 'study_name'
)

# 显示Top 10结果
print("\n===== Top 10 Predicted Drugs =====")
print(top_results.head(10))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from adjustText import adjust_text

def plot_drug_volcano(all_results, top_n_label=20,
                       prob_threshold=None,
                       rank_threshold=500,
                       study_name="BCBM"):
    """
    药物筛选火山图：
    X轴 = logit得分（预测强度）
    Y轴 = -log10(1-概率）（统计置信度）
    颜色 = Top/Middle/Low 候选
    """
    
    df = all_results.copy()
    
    # 计算绘图指标
    df['neg_log_p'] = -np.log10(1 - df['probability'] + 1e-10)
    
    # 分类
    if prob_threshold is None:
        prob_threshold = df.nsmallest(rank_threshold, 'rank')['probability'].min()
    
    df['category'] = 'Low priority'
    df.loc[df['rank'] <= rank_threshold, 'category'] = 'Candidate'
    df.loc[df['rank'] <= 50, 'category']  = 'Top candidate'
    
    color_map = {
        'Top candidate': '#E64B35',
        'Candidate':     '#F39B7F',
        'Low priority':  '#ADB5BD'
    }
    size_map  = {
        'Top candidate': 80,
        'Candidate':     40,
        'Low priority':  15
    }
    alpha_map = {
        'Top candidate': 0.95,
        'Candidate':     0.7,
        'Low priority':  0.3
    }
    
    fig, ax = plt.subplots(figsize=(12, 9))
    
    # 分层绘制（先低后高，保证Top在最上层）
    for cat in ['Low priority', 'Candidate', 'Top candidate']:
        sub = df[df['category'] == cat]
        ax.scatter(
            sub['logit_score'], sub['neg_log_p'],
            c     = color_map[cat],
            s     = size_map[cat],
            alpha = alpha_map[cat],
            label = f"{cat} (n={len(sub)})",
            zorder = {'Low priority':1,'Candidate':2,'Top candidate':3}[cat],
            edgecolors = 'white' if cat != 'Low priority' else 'none',
            linewidths = 0.5
        )
    
    # 添加参考线
    ax.axhline(y=-np.log10(1 - prob_threshold + 1e-10),
               color='grey', linestyle='--', linewidth=1, alpha=0.7)
    ax.axvline(x=0, color='grey', linestyle='--', linewidth=0.8, alpha=0.5)
    
    # 标注 Top N 化合物
    top_df = df.nsmallest(top_n_label, 'rank')
    texts  = []
    for _, row in top_df.iterrows():
        texts.append(ax.text(
            row['logit_score'], row['neg_log_p'],
            row['compound_id'],
            fontsize = 7.5,
            fontweight = 'bold' if row['rank'] <= 20 else 'normal',
            color = '#E64B35' if row['rank'] <= 10 else '#333333'
        ))
    
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='-', color='grey', lw=0.5),
                expand_points=(1.5, 1.5))
    
    # 美化
    ax.set_xlabel('Logit Score（预测强度）', fontsize=13)
    ax.set_ylabel('-log₁₀(1 - Probability)（置信度）', fontsize=13)
    ax.set_title(
        f'DrugReflector 药物筛选火山图\n{study_name}: 乳腺癌脑转移',
        fontsize=14, fontweight='bold', pad=15
    )
    
    ax.legend(loc='upper left', framealpha=0.9, fontsize=10)
    ax.spines[['top','right']].set_visible(False)
    
    # 统计标注
    n_top    = (df['category'] == 'Top candidate').sum()
    n_cand   = (df['category'] == 'Candidate').sum()
    n_total  = len(df)
    ax.text(0.98, 0.02,
            f"总化合物: {n_total}\n"
            f"Top候选 (Rank≤50): {n_top}\n"
            f"候选 (Rank≤{rank_threshold}): {n_cand}",
            transform = ax.transAxes,
            fontsize  = 9,
            ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5',
                      facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f'./figures/drug_volcano_{study_name}.pdf',
                dpi=300, bbox_inches='tight')
    plt.savefig(f'./figures/drug_volcano_{study_name}.png',
                dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"火山图已保存")
    print(f"Top 10 候选化合物:")
    print(df.nsmallest(10, 'rank')[
        ['compound_id','rank','probability','logit_score']
    ].to_string(index=False))
    
    return df

drug_volcano = plot_drug_volcano(
    all_results,
    top_n_label    = 20,
    rank_threshold = 500,
    study_name     = "BC"
)

plot_df = all_results.copy()

# Reset the index to make the compound ID a column
plot_df.reset_index(inplace=True)
plot_df.rename(columns={'index': 'compound_id'}, inplace=True)

# Rename columns to match what the plotting function expects
# Find the actual column names from your DataFrame
rank_col = [col for col in plot_df.columns if 'rank_' in col][0]
prob_col = [col for col in plot_df.columns if 'prob_' in col][0]

plot_df.rename(columns={
    rank_col: 'rank',
    'score': 'logit_score', # We named the logit score 'score' in the previous step
    prob_col: 'probability'
}, inplace=True)
# --- End of Data Preparation ---


# Now, you can run your plotting function without any changes to it.
# Make sure the 'plot_drug_volcano' function is defined in a cell above this one.
drug_volcano_df = plot_drug_volcano(
    all_results=plot_df, # Use the prepared DataFrame
    top_n_label=20,
    rank_threshold=500,
    study_name="BC"
)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from adjustText import adjust_text # This line is the important one
import os
def plot_drug_volcano(all_results, top_n_label=20,
                      prob_threshold=None,
                      rank_threshold=500,
                      study_name="BCBM"):
    """
    药物筛选火山图
    """
    
    # Make sure the output directory exists
    output_dir = "./figures"
    os.makedirs(output_dir, exist_ok=True)
    
    df = all_results.copy()
    
    # The rest of your function is perfect, no changes needed here
    # ... (all the plotting logic from your original code) ...
    
    df['neg_log_p'] = -np.log10(df['probability'] + 1e-10)
 
    if prob_threshold is None:
        prob_threshold = df.nsmallest(rank_threshold, 'rank')['probability'].min()
    
    df['category'] = 'Low priority'
    df.loc[df['rank'] <= rank_threshold, 'category'] = 'Candidate'
    df.loc[df['rank'] <= 50, 'category'] = 'Top candidate'
    
    color_map = {
        'Top candidate': '#E64B35', 'Candidate': '#F39B7F', 'Low priority': '#ADB5BD'
    }
    size_map = {
        'Top candidate': 80, 'Candidate': 40, 'Low priority': 15
    }
    alpha_map = {
        'Top candidate': 0.95, 'Candidate': 0.7, 'Low priority': 0.3
    }
    
    fig, ax = plt.subplots(figsize=(12, 9))
    
    for cat in ['Low priority', 'Candidate', 'Top candidate']:
        sub = df[df['category'] == cat]
        ax.scatter(
            sub['logit_score'], sub['neg_log_p'], c=color_map[cat], s=size_map[cat],
            alpha=alpha_map[cat], label=f"{cat} (n={len(sub)})",
            zorder={'Low priority':1,'Candidate':2,'Top candidate':3}[cat],
            edgecolors='white' if cat != 'Low priority' else 'none', linewidths=0.5
        )
    
    ax.axhline(y=-np.log10(1 - prob_threshold + 1e-10), color='grey', linestyle='--', linewidth=1, alpha=0.7)
    ax.axvline(x=0, color='grey', linestyle='--', linewidth=0.8, alpha=0.5)
    
    top_df = df.nsmallest(top_n_label, 'rank')
    texts = []
    for _, row in top_df.iterrows():
        texts.append(ax.text(
            row['logit_score'], row['neg_log_p'], row['compound_id'],
            fontsize=7.5, fontweight='bold' if row['rank'] <= 20 else 'normal',
            color='#E64B35' if row['rank'] <= 10 else '#333333'
        ))
    
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='-', color='grey', lw=0.5),
                expand_points=(1.5, 1.5))
    
    ax.set_xlabel('Logit Score (Prediction Strength)', fontsize=13)
    ax.set_ylabel('-log₁₀(1 - Probability) (Confidence)', fontsize=13)
    ax.set_title(f'DrugReflector Volcano Plot\n{study_name}: Breast Cancer Brain Metastasis',
                 fontsize=14, fontweight='bold', pad=15)
    
    ax.legend(loc='upper left', framealpha=0.9, fontsize=10)
    ax.spines[['top','right']].set_visible(False)
    
    n_top = (df['category'] == 'Top candidate').sum()
    n_cand = (df['category'] == 'Candidate').sum()
    n_total = len(df)
    ax.text(0.98, 0.02,
            f"Total Compounds: {n_total}\n"
            f"Top Candidates (Rank≤50): {n_top}\n"
            f"Candidates (Rank≤{rank_threshold}): {n_cand}",
            transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'drug_volcano_{study_name}.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, f'drug_volcano_{study_name}.png'), dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Volcano plot saved to {output_dir}/")
    print(f"Top 10 Candidate Compounds:")
    print(df.nsmallest(10, 'rank')[['compound_id','rank','probability','logit_score']].to_string(index=False))
    
    return df
drug_volcano_df = plot_drug_volcano(
    all_results=plot_df, # Use the prepared DataFrame
    top_n_label=20,
    rank_threshold=500,
    study_name="BC"
)