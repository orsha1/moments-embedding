import subprocess
import json
import os
import shutil
import sys


def create_dir(directory):
    os.makedirs(directory, exist_ok=True)


def copy_files(file_mappings):
    for source, destination in file_mappings:
        try:
            shutil.copyfile(source, destination)
        except FileNotFoundError:
            print("File not found:", source)
        except PermissionError:
            print("Permission denied for:", source)
        except Exception as e:
            print(f"An error occurred: {e}")


def main(config_file):
    main_dir = os.path.dirname(config_file)

    with open(config_file, "r") as file:
        config = json.load(file)

    if config["calculation"] == "step":
        config.update({"embed": "false", "quick": "true"})
        if not os.path.exists(os.path.join(main_dir, "embed")):
            print("Error: Embedding missing.")
            return
    else:
        config["embed"] = "true"
        config["quick"] = "true" if config["calculation"] == "full" else "false"
        create_dir(os.path.join(main_dir, "embed"))

    file_mappings = [(os.path.join(main_dir, config["in_lp"]), os.path.join(main_dir, "embed", "input_lattice_parameters.csv")), (os.path.join(main_dir, config["in_coords"]), os.path.join(main_dir, "embed", "input_coordinates.csv"))]
    copy_files(file_mappings)

    config["runpath"] = os.path.abspath(main_dir)
    config["datapath"] = os.path.join(config["runpath"], "embed")

    with open(os.path.join(main_dir, f"{config['name']}.in"), "w") as f:
        f.write(f"datapath={config['datapath']}/\n")
        f.write(f"runpath={config['runpath']}/\n")
        f.write(f"n_pbc={config['n_pbc']}\n")
        f.write(f"l_neis={config['l_neis']}\n")
        f.write(f"l_path={config['l_path']}\n")
        f.write(f"type={config['type']}\n")
        f.write(f"quick={config['quick']}\n")
        f.write(f"embed={config['embed']}")

    subprocess.run([os.path.join(config["embedpath"], "main"), os.path.join(main_dir, f"{config['name']}.in")])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <config_file>")
        sys.exit(1)
    config_file = sys.argv[1]
    if "/" not in config_file:
        config_file = f"./{config_file}"
    main(config_file)
