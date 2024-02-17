import csv 
from pathlib import Path
import re

import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(precision=15)
class Comparer:

    def __init__(self):
        self.read_folder_path = str(Path(__file__).parent) + "/results/"

    def run_prog(self):
        self.ask_number_of_files()
        self.ask_unit()
        self.initialize_graph()
        for i in range(self.number_of_files):
            while True:
                self.ask_name_of_files(i + 1)
                try:
                    with open(self.read_folder_path + self.file_name, "r") as file:
                        reader = csv.reader(file)

                        # Allocate memory
                        self.sol_time = np.zeros(50000)
                        self.energy = np.zeros(50000)

                        i = 0
                        for row in reader:
                            self.sol_time[i] = row[0]
                            self.energy[i] = row[1]
                            i += 1

                            # Extending memory buffer
                            if i % 50000 == 0:
                                self.sol_time = np.concatenate(self.sol_time, np.zeros(50000))
                                self.energy = np.concatenate(self.energy, np.zeros(50000))

                        self.sol_time = self.sol_time[:i]
                        self.energy = self.energy[:i]


                        self.ask_label_name()
                        self.plot_rel_energy()
                        break

                except FileNotFoundError:
                    print("File does not exit. Please try again.")

        self.show_plot()

    def ask_number_of_files(self):
        while True:
            try:
                self.number_of_files = int(input("Enter number of files you want to compare: ").strip())
                break
            except ValueError:
                print("Invalid input. Please try again.")

    def ask_unit(self):
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
        self.ax.set_title("Relative energy error against time")
        self.ax.set_xlabel(f"Time ({self.unit})")
        self.ax.set_ylabel("|(E(t)-E0)/E0|")

    def ask_name_of_files(self, file_number):
        while True:
            try:
                self.file_name = input(f"Enter the full name of file {file_number}(including .csv): ").strip()
                break
            except ValueError:
                print("Invalid input. Please try again.")

    def ask_label_name(self):
        self.label_name = input("Enter the label for this data: ").strip()

    def plot_rel_energy(self):
        self.ax.semilogy(
            self.sol_time,
            np.abs(
                (self.energy - self.energy[0])
                / self.energy[0]
            ),
            label=self.label_name,
        )

    def show_plot(self):
        self.fig.legend(loc=7)
        self.fig.tight_layout()
        self.fig.subplots_adjust(right=0.8)
        plt.show()


if __name__ == "__main__":
    comparer = Comparer()
    comparer.run_prog()