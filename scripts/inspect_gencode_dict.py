#!/usr/bin/env python3
"""
检查 Gencode 字典文件的内容和结构

用法:
    python scripts/inspect_gencode_dict.py /path/to/genecode_plus.dic
"""

import pickle
import sys
from pathlib import Path


def inspect_pickle_file(file_path: Path):
    """检查 pickle 文件的内容"""
    print(f"检查文件: {file_path}")
    print("=" * 80)

    try:
        with open(file_path, 'rb') as f:
            data = pickle.load(f)

        print(f"数据类型: {type(data)}")
        print(f"数据大小: {len(data) if hasattr(data, '__len__') else 'N/A'}")

        if isinstance(data, dict):
            print(f"\n字典键数量: {len(data)}")
            print(f"\n前5个键:")
            for i, key in enumerate(list(data.keys())[:5]):
                print(f"  {i+1}. {key}")
                value = data[key]
                print(f"     类型: {type(value)}")
                if isinstance(value, dict):
                    print(f"     子键: {list(value.keys())[:5]}")
                elif isinstance(value, (list, tuple)):
                    print(f"     长度: {len(value)}")
                    if len(value) > 0:
                        print(f"     第一个元素: {value[0]}")
                else:
                    print(f"     值: {str(value)[:100]}")
                print()

        elif isinstance(data, (list, tuple)):
            print(f"\n列表/元组长度: {len(data)}")
            print(f"前3个元素:")
            for i, item in enumerate(data[:3]):
                print(f"  {i+1}. {type(item)}: {str(item)[:100]}")

        else:
            print(f"\n数据内容: {str(data)[:500]}")

        print("\n" + "=" * 80)
        print("✓ 文件检查完成")

    except Exception as e:
        print(f"✗ 错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python inspect_gencode_dict.py <pickle_file>")
        sys.exit(1)

    file_path = Path(sys.argv[1])
    if not file_path.exists():
        print(f"错误: 文件不存在: {file_path}")
        sys.exit(1)

    inspect_pickle_file(file_path)
