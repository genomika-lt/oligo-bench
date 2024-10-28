import csv
import os
import re
import sys

from PyQt6.QtWidgets import (
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

        upper_table_buttons_layout = QHBoxLayout()
        save_yaml_button = QPushButton("Save YAML", self)
        save_yaml_button.clicked.connect(self.save_yaml)
        upper_table_buttons_layout.addWidget(save_yaml_button)
        layout.addLayout(upper_table_buttons_layout)

        # SECTION: Experiments (Table Widget)
        layout.addWidget(QLabel("Experiments"))
        self.experiments_table = QTableWidget(self)
        self.experiments_table.setColumnCount(2)  # path_to_sample, path_to_reference
        self.experiments_table.setHorizontalHeaderLabels(["Path to Sample", "Path to Reference"])
        self.experiments_table.horizontalHeader().setStretchLastSection(True)
        self.experiments_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.experiments_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)  # Disable default editing
        layout.addWidget(self.experiments_table)

        table_buttons_layout = QHBoxLayout()

        # Buttons below the table
        select_directory_button = QPushButton("Select Directory to Scan", self)
        select_directory_button.clicked.connect(self.select_scan_directory)
        table_buttons_layout.addWidget(select_directory_button)

        save_csv_button = QPushButton("Save Experiments CSV", self)
        save_csv_button.clicked.connect(self.save_experiments_csv)
        table_buttons_layout.addWidget(save_csv_button)

        layout.addLayout(table_buttons_layout)

        # Buttons for adding and deleting rows
        add_row_button = QPushButton("Add Row", self)
        add_row_button.clicked.connect(self.add_row)
        delete_row_button = QPushButton("Delete Row", self)
        delete_row_button.clicked.connect(self.delete_row)
        table_buttons_layout.addWidget(add_row_button)
        table_buttons_layout.addWidget(delete_row_button)

        self.experiments_table.cellDoubleClicked.connect(self.open_file_dialog)

        self.setLayout(layout)
        self.setWindowTitle('YAML and CSV File Generator')

    def select_samples_path(self):
        """
        Allows the user to select a CSV file for the experiment sheet.
        :return: None
        """
        try:
            file_path, _ = QFileDialog.getOpenFileName(self, "Select Experiment Sheet", "",
                                                       "CSV Files (*.csv);;All Files (*)")
            if file_path:
                self.samples_path.setText(file_path)
        except (OSError, IOError) as e:
            self.show_error(f"File I/O error selecting sample path: {str(e)}")
        except Exception as e:
            self.show_error(f"Unexpected error: {str(e)}")

    def select_output_directory(self):
        """
        Allows the user to select a directory for output files, including hidden folders.
        :return: None
        """
        try:
            options = QFileDialog.Option.DontUseNativeDialog | QFileDialog.Option.ShowDirsOnly
            folder_path = QFileDialog.getExistingDirectory(self,
                                                           "Select Output Directory",
                                                           options=options)
            if folder_path:
                self.output_folder.setText(folder_path)
        except (OSError, IOError) as e:
            self.show_error(f"File I/O error selecting output directory: {str(e)}")
        except Exception as e:
            self.show_error(f"Unexpected error: {str(e)}")

    def select_scan_directory(self):
        """
        Allows the user to select a directory to scan for experiment files, including hidden folders.
        :return: None
        """
        try:
            options = QFileDialog.Option.DontUseNativeDialog | QFileDialog.Option.ShowDirsOnly
            directory = QFileDialog.getExistingDirectory(self,
                                                         "Select Directory to Scan",
                                                         "", options=options)
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
        Shows an error if no matching folders are found.
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
            matches_found = False
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
                        matches_found = True

            if not matches_found:
                self.show_error("No matching experiment folders found in the selected directory.")
                return

            self.select_reference_file()

        except ValueError as e:
            self.show_error(f"Scan error: {str(e)}")
        except (OSError, IOError) as e:
            self.show_error(f"File system error during scan: {str(e)}")
        except Exception as e:
            self.show_error(f"Unexpected error during scan: {str(e)}")

    def select_reference_file(self):
        """
        Scans the selected directory for a reference file matching the pattern.
        If no valid .fa file is found, prompts the user to manually select it.
        :return: None
        """
        try:
            if not self.selected_directory:
                raise ValueError("No directory selected for scanning.")

            fa_pattern = re.compile(r"(test[/\\]resources[/\\]ref\.fa)$")
            found_fa_path = None

            for root, dirs, files in os.walk(self.selected_directory):
                for file in files:
                    if fa_pattern.search(os.path.join(root, file)):
                        potential_fa_path = os.path.join(root, file)
                        if os.path.getsize(potential_fa_path) > 0:
                            found_fa_path = potential_fa_path
                            break
                if found_fa_path:
                    break

            if found_fa_path:
                self.reference_file = found_fa_path
                for row in range(self.experiments_table.rowCount()):
                    self.experiments_table.setItem(row, 1, QTableWidgetItem(self.reference_file))
                QMessageBox.information(self, "Reference File Found",
                                        f"Reference file automatically found at {self.reference_file}.")

            else:
                QMessageBox.warning(self, "Reference File Not Found",
                                    "No valid reference file found in 'test/resources/'. "
                                    "Please select a .fa file manually.")
                file_path, _ = QFileDialog.getOpenFileName(self, "Select Reference File", "",
                                                           "FA Files (*.fa);;All Files (*)")
                if file_path:
                    if os.path.getsize(file_path) == 0:
                        self.show_error("The selected .fa file is empty. Please select a valid file.")
                        return
                    self.reference_file = file_path
                    for row in range(self.experiments_table.rowCount()):
                        self.experiments_table.setItem(row, 1, QTableWidgetItem(self.reference_file))

        except ValueError as e:
            self.show_error(f"Selection error: {str(e)}")
        except (OSError, IOError) as e:
            self.show_error(f"File system error during selection: {str(e)}")
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

            if not samples:
                self.show_error("Path to experiment sheet cannot be empty.")
                return
            if not output_folder:
                self.show_error("Output directory cannot be empty.")
                return
            if not report_name:
                self.show_error("Report name cannot be empty.")
                return

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

            file_path, _ = QFileDialog.getSaveFileName(self, "Save YAML File", "",
                                                       "YAML Files (*.yaml);;All Files (*)")
            if not file_path:
                self.show_error("No file path specified. YAML file not saved.")
                return

            if not file_path.endswith('.yaml'):
                file_path += '.yaml'

            with open(file_path, 'w', encoding='utf-8') as file:
                file.write(yaml_content)

            QMessageBox.information(self, "Success", "YAML file saved successfully.")

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
            if self.experiments_table.rowCount() == 0:
                self.show_error("No experiments to save.")
                return

            output_folder = self.output_folder.text()

            if not output_folder:
                self.show_error("Output directory cannot be empty.")
                return

            experiments_data = []
            for row in range(self.experiments_table.rowCount()):
                path_to_sample = self.experiments_table.item(row, 0)
                path_to_reference = self.experiments_table.item(row, 1)
                if path_to_sample and path_to_reference:
                    experiments_data.append([
                        path_to_sample.text(),
                        path_to_reference.text()
                    ])

            file_path, _ = QFileDialog.getSaveFileName(self, "Save Experiments CSV", "",
                                                       "CSV Files (*.csv);;All Files (*)")
            if file_path:
                if not file_path.endswith('.csv'):
                    file_path += '.csv'

                with open(file_path, 'w', newline='', encoding='utf-8') as file:
                    writer = csv.writer(file)
                    writer.writerow(["path_to_sample", "path_to_reference"])
                    writer.writerows(experiments_data)

                QMessageBox.information(self, "Success", "Experiments CSV file saved successfully.")

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

    def add_row(self):
        """
        Adds a new row to the experiments table with file selection buttons.
        :return: None
        """
        row_position = self.experiments_table.rowCount()
        self.experiments_table.insertRow(row_position)

    def delete_row(self):
        """
        Deletes the selected row from the experiments table.
        :return: None
        """
        selected_row = self.experiments_table.currentRow()
        if selected_row >= 0:
            self.experiments_table.removeRow(selected_row)
        else:
            self.show_error("Please select a row to delete.")

    def create_row_widget(self):
        """
        Creates a custom widget for each row with file path input and button.
        :return: QWidget containing path input and button
        """
        row_widget = QWidget()
        layout = QHBoxLayout()
        path_input = QLineEdit(self)
        layout.addWidget(path_input)
        row_widget.setLayout(layout)
        return row_widget

    def open_file_dialog(self, row, column):
        """
        Opens a file dialog to select a file or directory and updates the specified cell.
        :param row: The row index of the cell that was double-clicked
        :param column: The column index of the cell that was double-clicked
        :return: None
        """
        if column == 0:
            folder_path = QFileDialog.getExistingDirectory(self, "Select Sample Folder", "")
            if folder_path:
                self.experiments_table.setItem(row, column, QTableWidgetItem(folder_path))
        elif column == 1:
            fa_file_path, _ = QFileDialog.getOpenFileName(self, "Select Reference File", "",
                                                          "FA Files (*.fa);;All Files (*)")
            if fa_file_path:
                self.experiments_table.setItem(row, column, QTableWidgetItem(fa_file_path))

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = YamlForm()
    ex.show()
    sys.exit(app.exec())
