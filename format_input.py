#!/usr/bin/env python3
"""
ASTRA / Generator 输入文件格式规范化工具。

功能：
  - CRLF → LF 换行符转换
  - 统一缩进（2 空格）
  - 去除行尾空白和多余空行
  - 统一 namelist 格式（key=value, 后加逗号）
  - 保留注释行和 `...` 分隔符

用法：
  python format_input.py <file>              # 原地格式化
  python format_input.py <file> --check      # 仅检查不修改
  python format_input.py <file> -o out.in    # 输出到指定文件
  python format_input.py --all               # 批量格式化所有 .in 文件
"""

import re
import sys
import argparse
from pathlib import Path
from typing import Optional


def normalize_line_endings(text: str) -> str:
    """CRLF → LF"""
    return text.replace("\r\n", "\n").replace("\r", "\n")


def strip_trailing_spaces(line: str) -> str:
    """去除行尾空白"""
    return line.rstrip()


def format_namelist_line(line: str, in_namelist: bool) -> str:
    """
    格式化 namelist 中的参数行。
    规则：
      - 统一 `key = value` → `key=value`
      - 如果行尾没有逗号且不是注释/块标记，添加逗号
    """
    stripped = line.strip()
    if not stripped:
        return ""

    # 保留注释行
    if stripped.startswith("!"):
        return "  " + stripped

    # 保留块标记 &NAME 和 /
    if stripped.startswith("&") or stripped == "/":
        return stripped

    # 格式化参数行
    # 去除 key 和 = 之间的空格
    formatted = re.sub(r"\s*=\s*", "=", stripped)
    # 去除参数名后 ( 前的空格（仅非引号内的）：WK_Z  (1) → WK_Z(1)
    if "'" not in formatted and '"' not in formatted:
        formatted = re.sub(r"\s+\(", "(", formatted)

    # 如果行尾没有逗号，添加
    if not formatted.endswith(","):
        formatted += ","

    return "  " + formatted


def format_input_file(content: str) -> str:
    """
    格式化整个输入文件内容。

    Args:
        content: 原始文件内容

    Returns:
        格式化后的内容
    """
    content = normalize_line_endings(content)
    lines = content.split("\n")

    result = []
    in_namelist = False
    prev_blank = False

    for line in lines:
        stripped = strip_trailing_spaces(line)

        # 检测 namelist 开始
        if stripped.strip().startswith("&"):
            in_namelist = True
            # 块前空行
            if result and result[-1] != "":
                result.append("")
            result.append(stripped)
            prev_blank = False
            continue

        # 检测 namelist 结束
        if in_namelist and stripped.strip() == "/":
            result.append(" /")
            result.append("")
            in_namelist = False
            prev_blank = True
            continue

        # namelist 内部
        if in_namelist:
            formatted = format_namelist_line(stripped, in_namelist)
            if formatted:
                result.append(formatted)
                prev_blank = False
            continue

        # namelist 外部（注释、... 分隔符等）
        if stripped:
            if stripped == "...":
                if result and result[-1] != "":
                    result.append("")
                result.append(stripped)
                result.append("")
                prev_blank = True
            else:
                result.append(stripped)
                prev_blank = False
        else:
            if not prev_blank and result:
                result.append("")
                prev_blank = True

    # 去除末尾多余空行
    while result and result[-1] == "":
        result.pop()

    return "\n".join(result) + "\n"


def format_file(filepath: Path, output: Optional[Path] = None, check_only: bool = False) -> bool:
    """
    格式化单个文件。

    Returns:
        True 如果文件已经被规范化（无需修改），False 如果需要修改
    """
    try:
        original = filepath.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        print(f"  ⊘ {filepath.name} — 跳过（非文本文件）")
        return True
    formatted = format_input_file(original)

    if original == formatted:
        print(f"  ✓ {filepath.name} — 已规范")
        return True

    if check_only:
        print(f"  ✗ {filepath.name} — 需要格式化")
        return False

    target = output or filepath
    target.write_text(formatted, encoding="utf-8")
    print(f"  ✓ {filepath.name} — 已格式化")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="ASTRA / Generator 输入文件格式规范化工具"
    )
    parser.add_argument(
        "file", nargs="?", help="要格式化的 .in 文件路径"
    )
    parser.add_argument(
        "-o", "--output", help="输出文件路径（默认原地修改）"
    )
    parser.add_argument(
        "--check", action="store_true", help="仅检查，不修改文件"
    )
    parser.add_argument(
        "--all", action="store_true", help="批量格式化项目中的所有 .in 文件"
    )
    args = parser.parse_args()

    if args.all:
        project_dir = Path.cwd()
        in_files = list(project_dir.rglob("*.in"))
        print(f"找到 {len(in_files)} 个 .in 文件\n")
        ok = 0
        for f in sorted(in_files):
            if format_file(f, check_only=args.check):
                ok += 1
        print(f"\n规范: {ok}/{len(in_files)}")
        return

    if not args.file:
        parser.print_help()
        return

    filepath = Path(args.file)
    if not filepath.exists():
        print(f"错误: 文件不存在 — {filepath}")
        sys.exit(1)

    output = Path(args.output) if args.output else None
    format_file(filepath, output=output, check_only=args.check)


if __name__ == "__main__":
    main()
