'''将 centroids_all_alphafold.csv 和 centroids_all_pdb.csv 两个文件的内容整合为新的 all_pockets.txt。输出格式与附件一致，
并且会保留 alphafold 文件中的蛋白质名称（例如 A0A1B1WA13_alphafold）以及 pdb 文件中的原名（例如 F9UM18_7EBO）。脚本与这三个文件在同一目录下运行
成功'''
import csv
import os

def format_center(value: str) -> str:
    """Format numeric string to one decimal place, matching the example output."""
    try:
        v = float(value)
    except ValueError:
        # If parsing fails, keep original text (unlikely given CSV format)
        return value
    # Avoid "-0.0"
    v = 0.0 if abs(v) < 5e-16 else v
    return f"{v:.1f}"

def process_csv(input_csv: str, out_lines: list) -> None:
    """Read a CSV file and append formatted lines to out_lines."""
    with open(input_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["source_pdb"]  # Use as-is for both alphafold and pdb
            x = format_center(row["centroid_x"])
            y = format_center(row["centroid_y"])
            z = format_center(row["centroid_z"])
            line = (
                f"{name}\t"
                f"--center_x {x} --center_y {y} --center_z {z} "
                f"--size_x 20.0 --size_y 20.0 --size_z 20.0"
            )
            out_lines.append(line)

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    input_files = ["centroids_all_alphafold.csv", "centroids_all_pdb.csv"]
    lines = []
    for filename in input_files:
        path = os.path.join(script_dir, filename)
        process_csv(path, lines)

    # Write all lines to all_pockets.txt (in the script directory)
    out_path = os.path.join(script_dir, "all_pockets.txt")
    with open(out_path, "w", encoding="utf-8", newline="\n") as out:
        for line in lines:
            out.write(line + "\n")

if __name__ == "__main__":
    main()