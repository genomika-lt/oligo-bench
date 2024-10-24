import sys
import csv
import os
import re
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QPushButton, QFileDialog, QCheckBox, QTableWidget,
    QTableWidgetItem, QHeaderView, QAbstractItemView, QMessageBox
)

class YamlForm(QWidget):
    def __init__(self):
        """
        Initializes the main form and UI elements.
        """
        super().__init__()

        self.selected_directory = ""
        self.reference_file = ""
        self.initUI()

    def initUI(self):
        """
        Sets up the user interface for the form.
        :return: None
        """
        layout = QVBoxLayout()

        # SECTION: Samples (Button and Text Field)
        samples_layout = QHBoxLayout()
        self.samples_path = QLineEdit(self)
        samples_layout.addWidget(self.samples_path)
        samples_button = QPushButton("Select Experiment Sheet", self)
        samples_button.clicked.connect(self.select_samples_path)
        samples_layout.addWidget(samples_button)
        layout.addWidget(QLabel("Path to experiment sheet (CSV format):"))
        layout.addLayout(samples_layout)

        # SECTION: Output (Button and Text Field)
        output_layout = QHBoxLayout()
        self.output_folder = QLineEdit(self)
        output_layout.addWidget(self.output_folder)
        output_button = QPushButton("Select Output Directory", self)
        output_button.clicked.connect(self.select_output_directory)
        output_layout.addWidget(output_button)
        layout.addWidget(QLabel("Output directory:"))
        layout.addLayout(output_layout)

        # SECTION: Report Name
        self.report_name = QLineEdit(self)
        layout.addWidget(QLabel("Report Name:"))
        layout.addWidget(self.report_name)

        # SECTION: Basecalling (Checkbox and Dropdown)
        self.basecall_checkbox = QCheckBox("Enable Basecall", self)
        layout.addWidget(self.basecall_checkbox)
        layout.addWidget(QLabel("Dorado Model:"))
        self.dorado_model_combobox = QComboBox(self)
        self.dorado_model_combobox.addItems(["fast", "hac", "sup"])
        layout.addWidget(self.dorado_model_combobox)

        # SECTION: Experiments (Table Widget)
        layout.addWidget(QLabel("Experiments (path_to_sample, path_to_reference):"))

        self.experiments_table = QTableWidget(self)
        self.experiments_table.setColumnCount(2)  # path_to_sample, path_to_reference
        self.experiments_table.setHorizontalHeaderLabels(["Path to Sample", "Path to Reference"])
        self.experiments_table.horizontalHeader().setStretchLastSection(True)
        self.experiments_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.experiments_table.setEditTriggers(QAbstractItemView.DoubleClicked)
        layout.addWidget(self.experiments_table)

        # Buttons below the table
        table_buttons_layout = QHBoxLayout()
        select_directory_button = QPushButton("Select Directory to Scan", self)
        select_directory_button.clicked.connect(self.select_scan_directory)
        save_yaml_button = QPushButton("Save YAML", self)
        save_yaml_button.clicked.connect(self.save_yaml)
        save_csv_button = QPushButton("Save Experiments CSV", self)
        save_csv_button.clicked.connect(self.save_experiments_csv)
        table_buttons_layout.addWidget(select_directory_button)
        table_buttons_layout.addWidget(save_yaml_button)
        table_buttons_layout.addWidget(save_csv_button)
        layout.addLayout(table_buttons_layout)

        self.setLayout(layout)
        self.setWindowTitle('YAML and CSV File Generator')

    def select_samples_path(self):
        """
        Allows the user to select a CSV file for the experiment sheet.
        :return: None
        """
        try:
            options = QFileDialog.Options()
            file_path, _ = QFileDialog.getOpenFileName(self, "Select Experiment Sheet", "",
                                                       "CSV Files (*.csv);;All Files (*)",
                                                       options=options)
            if file_path:
                self.samples_path.setText(file_path)
        except (OSError, IOError) as e:
            self.show_error(f"File I/O error selecting sample path: {str(e)}")
        except Exception as e:
            self.show_error(f"Unexpected error: {str(e)}")

    def select_output_directory(self):
        """
        Allows the user to select a directory for output files.
        :return: None
        """
        try:
            options = QFileDialog.Options()
            folder_path = QFileDialog.getExistingDirectory(self, "Select Output Directory", options=options)
            if folder_path:
                self.output_folder.setText(folder_path)
        except (OSError, IOError) as e:
            self.show_error(f"File I/O error selecting output directory: {str(e)}")
        except Exception as e:
            self.show_error(f"Unexpected error: {str(e)}")

    def select_scan_directory(self):
        """
        Allows the user to select a directory to scan for experiment files.
        :return: None
        """
        try:
            options = QFileDialog.Options()
            directory = QFileDialog.getExistingDirectory(self, "Select Directory to Scan", options=options)
            if directory:
                self.selected_directory = directory
                self.scan_experiment_files()
        except (OSError, IOError) as e:
            self.show_error(f"Directory selection error: {str(e)}")
        except Exception as e:
            self.show_error(f"Unexpected error: {str(e)}")

    def scan_experiment_files(self):
        """
        Scans the selected directory for experiment files matching a specific pattern.
        :return: None
        """
        try:
            if not self.selected_directory:
                raise ValueError("No directory selected for scanning.")

            pattern = re.compile(
                r'(test_experiment_\d+[\\/]+test_sample_\d+[\\/]+20\d{2}(0[1-9]|'
                r'1[0-2])(0[1-9]|[12]\d|3[01])_\d{4}_MN\d{5}_ATQ\d+_\w+)$'

            )

            self.experiments_table.setRowCount(0)

            for root, dirs, files in os.walk(self.selected_directory):
                for dir_name in dirs:
                    full_path = os.path.join(root, dir_name)
                    match = pattern.search(full_path)
                    if match:
                        sample_path = os.path.join(root, match.group(0))
                        row_position = self.experiments_table.rowCount()
                        self.experiments_table.insertRow(row_position)
                        self.experiments_table.setItem(row_position, 0, QTableWidgetItem(sample_path))
                        self.experiments_table.setItem(row_position, 1, QTableWidgetItem(""))

            self.select_reference_file()
        except ValueError as e:
            self.show_error(f"Scan error: {str(e)}")
        except (OSError, IOError) as e:
            self.show_error(f"File system error during scan: {str(e)}")
        except Exception as e:
            self.show_error(f"Unexpected error during scan: {str(e)}")

    def select_reference_file(self):
        """
        Allows the user to select a reference file for the experiments.
        :return: None
        """
        try:
            options = QFileDialog.Options()
            file_path, _ = QFileDialog.getOpenFileName(self, "Select Reference File", "",
                                                       "FA Files (*.fa);;All Files (*)",
                                                       options=options)
            if file_path:
                self.reference_file = file_path
                for row in range(self.experiments_table.rowCount()):
                    self.experiments_table.item(row, 1).setText(self.reference_file)
        except (OSError, IOError) as e:
            self.show_error(f"File selection error: {str(e)}")
        except Exception as e:
            self.show_error(f"Unexpected error: {str(e)}")

    def save_yaml(self):
        """
        Saves the user input to a YAML file.
        :return: None
        """
        try:
            samples = self.samples_path.text()
            output_folder = self.output_folder.text()
            report_name = self.report_name.text()
            basecall = self.basecall_checkbox.isChecked()
            dorado_model = self.dorado_model_combobox.currentText()

            yaml_content = (
                "# DEFAULT\n"
                "# path to experiment sheet (CSV format)\n"
                f"samples: {samples}\n\n\n"
                "# OUTPUT\n"
                "# directory\n"
                f"output_folder: {output_folder}\n"
                "# Report Name\n"
                f"report_name: {report_name}\n\n\n"
                "# BASECALLING\n"
                "# define parameters for a specific rule\n"
                f"basecall: {'True' if basecall else 'False'}\n"
                "# dorado model 'fast', 'hac' or 'sup'\n"
                f"dorado_model: {dorado_model}\n"
            )

            options = QFileDialog.Options()
            file_path, _ = QFileDialog.getSaveFileName(self, "Save YAML File", "",
                                                       "YAML Files (*.yaml);;All Files (*)",
                                                       options=options)
            if file_path:
                with open(file_path, 'w',encoding='utf-8') as file:
                    file.write(yaml_content)
        except (OSError, IOError) as e:
            self.show_error(f"Error writing YAML file: {str(e)}")
        except Exception as e:
            self.show_error(f"Unexpected error saving YAML: {str(e)}")

    def save_experiments_csv(self):
        """
        Saves the experiments data to a CSV file.
        :return: None
        """
        try:
            experiments_data = []
            for row in range(self.experiments_table.rowCount()):
                path_to_sample = self.experiments_table.item(row, 0)
                path_to_reference = self.experiments_table.item(row, 1)
                if path_to_sample and path_to_reference:
                    experiments_data.append([
                        path_to_sample.text(),
                        path_to_reference.text()
                    ])

            options = QFileDialog.Options()
            file_path, _ = QFileDialog.getSaveFileName(self, "Save Experiments CSV", "",
                                                       "CSV Files (*.csv);;All Files (*)",
                                                       options=options)
            if file_path:
                with open(file_path, 'w', newline='',encoding='utf-8') as file:
                    writer = csv.writer(file)
                    writer.writerow(["path_to_sample", "path_to_reference"])
                    writer.writerows(experiments_data)
        except (OSError, IOError) as e:
            self.show_error(f"Error writing CSV file: {str(e)}")
        except Exception as e:
            self.show_error(f"Unexpected error saving CSV: {str(e)}")

    def show_error(self, message):
        """
        Displays an error message in a popup window.
        :param message: message to display
        :return: None
        """
        QMessageBox.critical(self, "Error", message)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = YamlForm()
    ex.show()
    sys.exit(app.exec_())
