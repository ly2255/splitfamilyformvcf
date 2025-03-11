import vcf
import argparse
from collections import defaultdict
import os
import chardet


def parse_pedigree(ped_file):
    """解析家系文件，返回每个家系的样本列表"""
    families = defaultdict(list)

    # 检测文件编码
    with open(ped_file, 'rb') as f:
        raw_data = f.read()
        result = chardet.detect(raw_data)
        encoding = result['encoding']
        print(f"Detected encoding: {encoding}")

    with open(ped_file, 'r', encoding=encoding) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split()
            if len(fields) >= 6:
                family_id, sample_id, father_id, mother_id, _, _ = fields[:6]

                # 添加proband
                families[family_id].append(sample_id)
                # 添加父母（排除0表示的未知父母）
                families[family_id].append(father_id)
                families[family_id].append(mother_id)

    return families


def split_vcf_by_family(vcf_file, families, outdir):
    """将VCF按家系分割为单独的VCF文件"""
    with open(vcf_file, 'r') as vcf_f:
        lines = vcf_f.readlines()

    # 提取VCF头部和 #CHROM 行
    header_lines = []
    chrom_line = None
    for line in lines:
        if line.startswith("##"):
            header_lines.append(line.strip())
        elif line.startswith("#CHROM"):
            chrom_line = line.strip().split("\t")
            break

    if chrom_line is None:
        raise ValueError("VCF 文件缺少 #CHROM 行！")

    sample_names = chrom_line[9:]  # 获取所有样本名

    # 确保输出目录存在
    os.makedirs(outdir, exist_ok=True)

    for family_id, members in families.items():
        # 过滤家系中的有效样本
        family_samples = [s for s in members if s in sample_names]

        if len(family_samples) < 1:
            print(f"❌ Skipping family {family_id}: Not enough valid samples in VCF")
            continue

        # 输出家系样本信息
        proband = members[0] if members[0] in sample_names else "N/A"
        father = members[1] if len(members) > 1 and members[1] in sample_names else "N/A"
        mother = members[2] if len(members) > 2 and members[2] in sample_names else "N/A"
        print(f"✅ Processing family {family_id}: Proband={proband}, Father={father}, Mother={mother}")

        output_file = os.path.join(outdir, f"{family_id}.vcf")

        with open(output_file, 'w') as out_f:
            # 写入头部注释行
            for header in header_lines:
                out_f.write(header + "\n")

            # 生成家系的 `#CHROM` 行
            selected_cols = chrom_line[:9] + family_samples
            out_f.write("\t".join(selected_cols) + "\n")

            # 逐行写入VCF变异信息
            for line in lines:
                if not line.startswith("#"):
                    fields = line.strip().split("\t")
                    genotypes = fields[9:]

                    # 只保留家系样本的数据
                    sample_data = {sample: genotypes[idx] for idx, sample in enumerate(sample_names)}
                    family_data = [sample_data[s] for s in family_samples]

                    # 重新组合行
                    output_line = fields[:9] + family_data
                    out_f.write("\t".join(output_line) + "\n")

        print(f"✅ Successfully wrote VCF for family {family_id}: {output_file}")

    print("✅ VCF splitting completed!")


def main():
    parser = argparse.ArgumentParser(description="Split VCF by family")
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--ped', required=True, help='Pedigree file')
    parser.add_argument('--outdir', required=True, help='Output directory')
    args = parser.parse_args()

    # 解析家系信息
    families = parse_pedigree(args.ped)
    print(f"📂 Loaded {len(families)} families from PED file")

    # 拆分 VCF
    print("📂 Splitting VCF by family...")
    split_vcf_by_family(args.vcf, families, args.outdir)

    print("✅ All families processed successfully!")


if __name__ == "__main__":
    main()
# python extractfamilyfromvcf.py --vcf all.sample.remaining.vcf --ped ped.txt --outdir ./family
# pedfile
#FamilyID	SampleID	FatherID	MotherID	Sex	Phenotype
#DY-ZYDF-300-052	DY-ZYDF-300-052	DY-ZYDF-298-063	DY-ZYDF-296-088	2	2
#DY-ZYDF-298-065	DY-ZYDF-298-065	DY-ZYDF-298-066	DY-ZYDF-298-067	2	2
#DY-ZYDF-298-079	DY-ZYDF-298-079	DY-ZYDF-298-080	DY-ZYDF-298-081	2	2
#DY-ZYDF-298-074	DY-ZYDF-298-074	DY-ZYDF-298-075	DY-ZYDF-298-076	2	2
#DY-ZYDF-296-098	DY-ZYDF-296-098	DY-ZYDF-296-099	DY-ZYDF-300-053	2	2
#DY-ZYDF-296-089	DY-ZYDF-296-089	DY-ZYDF-296-090	DY-ZYDF-296-091	2	2
#DY-ZYDF-298-068	DY-ZYDF-298-068	DY-ZYDF-298-069	DY-ZYDF-298-070	2	2
#DY-ZYDF-296-095	DY-ZYDF-296-095	DY-ZYDF-296-096	DY-ZYDF-296-097	2	2
#DY-ZYDF-300-054	DY-ZYDF-300-054	DY-ZYDF-298-077	DY-ZYDF-298-078	2	2
#DY-ZYDF-300-051	DY-ZYDF-300-051	DY-ZYDF-296-086	DY-ZYDF-296-087	2	2
#DY-ZYDF-296-092	DY-ZYDF-296-092	DY-ZYDF-298-064	DY-ZYDF-296-093	2	2
#DY-ZYDF-298-071	DY-ZYDF-298-071	DY-ZYDF-298-072	DY-ZYDF-298-073	2	2
#DY-ZYDF-300-050	DY-ZYDF-300-050	0	0	2	2
#DY-ZYDF-296-094	DY-ZYDF-296-094	0	0	2	2
