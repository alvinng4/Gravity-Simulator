"""
Compare the energy of the data inside the gravity_plot/results file

IMPORTANT: matplotlib cannot be interrupted by input() 
"""
import argparse
import csv
from pathlib import Path
import re
import sys

import numpy as np
import matplotlib.pyplot as plt


class Comparer:
    def __init__(self):
        self.read_folder_path = Path(__file__).parent / "results"
        self.read_command_line_arg()
        self.title = self.args.title
        self.file_names_labels = {}

    def run_prog(self):
        self.ask_unit()
        self.ask_files_info()
        self.initialize_graph()
        for file_name in self.file_names_labels:
            while True:
                try:
                    with open(self.read_folder_path / file_name, "r") as file:
                        reader = csv.reader(file)

                        # Allocate memory
                        self.sol_time = np.zeros(50000)
                        self.energy = np.zeros(50000)

                        i = 0
                        for row in reader:
                            self.sol_time[i] = row[0]
                            self.energy[i] = row[2]
                            i += 1

                            # Extending memory buffer
                            if i % 50000 == 0:
                                self.sol_time = np.concatenate(
                                    (self.sol_time, np.zeros(50000))
                                )
                                self.energy = np.concatenate(
                                    (self.energy, np.zeros(50000))
                                )

                        self.sol_time = self.sol_time[:i]
                        self.energy = self.energy[:i]

                        self.plot_rel_energy(self.file_names_labels[file_name])
                        break

                except FileNotFoundError:
                    sys.exit("Error: file is not found. Exiting the program")

        self.show_plot()

    def ask_unit(self):
        """
        Ask user the time unit (day/year) of all the data.
        """
        while True:
            self.unit = input("Enter the unit of time (day/year): ")
            if matches := re.search(
                r"^\s*(day|year|days|years|d|y)\s*$", self.unit, re.IGNORECASE
            ):
                if matches.group(1):
                    if matches.group(1).lower() in ["day", "days", "d"]:
                        self.unit = "days"
                        break
                    elif matches.group(1).lower() in ["year", "years", "y"]:
                        self.unit = "years"
                        break

            print("Invalid input. Please try again.")

    def initialize_graph(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title(self.title)
        self.ax.set_xlabel(f"Time ({self.unit})")
        self.ax.set_ylabel("|(E(t)-E0)/E0|")

    def ask_files_info(self):
        """
        Ask user the number of files, all the file names and their respective labels
        """
        while True:
            try:
                self.number_of_files = int(
                    input("Enter number of files you want to compare: ").strip()
                )
                break
            except ValueError:
                print("Invalid input. Please try again.")

        for i in range(self.number_of_files):
            while True:
                self.file_name = input(
                    f"Enter the full name of file {i + 1}(including .csv): "
                ).strip()
                self.file_path = self.read_folder_path / self.file_name
                if self.file_path.is_file():
                    self.label_name = input("Enter the label for this data: ").strip()
                    self.file_names_labels[self.file_name] = self.label_name
                    break
                else:
                    print("File does not exist! Please try again.")

    def plot_rel_energy(self, label_name):
        self.ax.semilogy(
            self.sol_time,
            np.abs((self.energy - self.energy[0]) / self.energy[0]),
            label=label_name,
        )

    def show_plot(self):
        self.fig.legend(loc=7)
        self.fig.tight_layout()
        self.fig.subplots_adjust(right=0.8)
        plt.show()

    def read_command_line_arg(self):
        parser = argparse.ArgumentParser(description="N-body gravity simulator")
        parser.add_argument(
            "--title",
            "-t",
            default="Relative energy error against time",
            type=str,
            help="Usage: --title <title>",
        )
        self.args = parser.parse_args()


if __name__ == "__main__":
    comparer = Comparer()
    comparer.run_prog()
