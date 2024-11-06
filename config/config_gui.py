"""
YamlForm Application

This module defines a PyQt6 application that provides
a graphical user interface (GUI) for managing YAML and CSV files
related to experiments. The application allows users
to configure base calling options, select dorado models, and
manage experiment paths through a table interface.
"""
import csv
import os
import sys
import yaml
from PyQt6.QtCore import Qt,QTimer # pylint: disable=E0611
from PyQt6.QtWidgets import (QApplication, QWidget, QVBoxLayout,
                             QHBoxLayout, QLabel,QComboBox, QPushButton,
                             QCheckBox, QTableWidget, QTableWidgetItem,
                             QHeaderView, QAbstractItemView, QMessageBox,
                             QTextEdit, QSplitter, QFileDialog)


class YamlForm(QWidget):
    """
    A QWidget subclass that provides a user interface for loading and managing
    YAML configuration files.

    This class allows users to interact with YAML configuration data, enabling
    them to load, display, modify, and save configurations via the UI. It provides
    various widgets for editing different aspects of the YAML file and ensures
    that changes can be saved back to the file.
    """
    def __init__(self):
        """
        Initializes the main form and UI elements.
        """
        super().__init__()
        self.elapsed_time = None
        self.unsaved_changes = None
        self.timer = None
        if getattr(sys, 'frozen', False) or "__compiled__" in globals():
            self.project_root = os.path.abspath(os.path.join(os.getcwd()))
        else:
            self.project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        self.initUI()
        self.unset_unsaved_changes()
        self.stop_requested = False
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.Window)
        self.showMaximized()

        yaml_path = os.path.join(self.project_root, "config", "config.yaml")
        csv_path = os.path.join(self.project_root, "config", "experiments.csv")

        if not self.check_permissions(yaml_path) or not self.check_permissions(csv_path):
            return

        self.load_csv_on_startup(csv_path)
        self.load_yaml_on_startup(yaml_path)
        self.running = False

    def initUI(self):
        """
        Sets up the user interface for the form.
        :return: None
        """
        # pylint: disable=C0103
        layout = QVBoxLayout()

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setHandleWidth(10)
        # SECTION: Left side - Controls and Log Window
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)

        # SECTION: Controls (Run button, Basecall checkbox, Dorado model combobox, Update button)
        update_button = QPushButton("Update", self)
        left_layout.addWidget(update_button)

        self.basecall_checkbox = QCheckBox("Enable Basecall", self)
        self.basecall_checkbox.stateChanged.connect(self.set_unsaved_changes)
        left_layout.addWidget(self.basecall_checkbox)

        left_layout.addWidget(QLabel("Dorado Model:"))
        self.dorado_model_combobox = QComboBox(self)
        self.dorado_model_combobox.addItems(["fast", "hac", "sup"])
        self.dorado_model_combobox.currentIndexChanged.connect(self.set_unsaved_changes)
        left_layout.addWidget(self.dorado_model_combobox)

        self.run_stop_button = QPushButton("Run", self)
        self.run_stop_button.clicked.connect(self.toggle_run_stop)
        left_layout.addWidget(self.run_stop_button)

        # SECTION: Log Window
        self.log_window = QTextEdit(self)
        self.log_window.setReadOnly(True)
        left_layout.addWidget(self.log_window)

        splitter.addWidget(left_widget)

        # SECTION: Right side - Experiments Table
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)

        right_layout.addWidget(QLabel("Experiments"))
        self.experiments_table = QTableWidget(self)
        self.experiments_table.setColumnCount(2)  # Path to Sample and Path to Reference
        self.experiments_table.setHorizontalHeaderLabels(["Path to Sample", "Path to Reference"])
        self.experiments_table.horizontalHeader().setStretchLastSection(True)
        (self.experiments_table.horizontalHeader().
         setSectionResizeMode(QHeaderView.ResizeMode.Stretch))
        self.experiments_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.experiments_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        right_layout.addWidget(self.experiments_table)

        # Buttons for adding and deleting rows
        table_buttons_layout = QHBoxLayout()
        add_row_button = QPushButton("Add", self)
        add_row_button.clicked.connect(self.add_row)
        table_buttons_layout.addWidget(add_row_button)

        delete_row_button = QPushButton("Delete", self)
        delete_row_button.clicked.connect(self.delete_row)
        table_buttons_layout.addWidget(delete_row_button)

        right_layout.addLayout(table_buttons_layout)

        splitter.addWidget(right_widget)

        layout.addWidget(splitter)

        splitter.setSizes([self.width() // 3, self.width() * 2 // 3])
        self.experiments_table.cellDoubleClicked.connect(self.open_file_dialog)

        self.setLayout(layout)
        self.setWindowTitle('YAML and CSV File Generator')
        self.load_stylesheet()

    def toggle_run_stop(self):
        """
        Toggles the Run/Stop button text and state,
         prompting to save CSV and YAML files if necessary.
        """
        if self.running:
            self.stop_requested = True
            self.run_stop_button.setText("Stopping...")
        else:
            csv_reply = QMessageBox.question(
                self,
                'Save CSV File',
                "Do you want to save the CSV file before running?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No
            )
            if csv_reply == QMessageBox.StandardButton.Yes:
                self.save_experiments_csv()

            yaml_reply = QMessageBox.question(
                self,
                'Save YAML File',
                "Do you want to save the YAML file before running?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No
            )
            if yaml_reply == QMessageBox.StandardButton.Yes:
                self.save_yaml()

            self.log_window.clear()
            self.run_stop_button.setText("Stop")
            self.run_stop_button.setStyleSheet("background-color: red; color: white;")
            self.unset_unsaved_changes()
            self.running = True
            self.stop_requested = False

            self.start_run()

    def start_run(self):
        """
        Simulates a long-running process for 10 seconds.
        """
        self.log_window.clear()
        self.elapsed_time = 0
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_log)
        self.timer.start(1000)

    def update_log(self):
        """
        Updates the log window every second and checks if the stop is requested.
        """
        if self.stop_requested:
            self.handle_stop()
            return
        self.elapsed_time += 1
        self.log_window.append(f"Running... {self.elapsed_time} second(s)")
        if self.elapsed_time >= 10:
            self.timer.stop()
            self.finish_run()

    def handle_stop(self):
        """
        Handles stopping the execution and updates the UI.
        """
        self.timer.stop()
        self.log_window.append("<font color='red'>Execution was stopped.</font>")
        self.run_stop_button.setText("Run")
        self.run_stop_button.setStyleSheet("")
        self.running = False
        self.stop_requested = False

    def finish_run(self):
        """
        Resets the state after the simulated run finishes.
        """
        self.timer.stop()
        self.running = False
        self.run_stop_button.setText("Run")
        self.run_stop_button.setStyleSheet("")
        self.log_window.append("Run finished.")

    def check_permissions(self, path):
        """
        Check if the application has read, write, and execute permissions for the specified path.
        :param path: Path to check
        :return: True if all permissions are granted, otherwise False
        """
        if not os.path.exists(path):
            self.show_error(f"File not found: {path}")
            return False

        permissions = {
            "read": os.access(path, os.R_OK),
            "write": os.access(path, os.W_OK),
            "execute": os.access(path, os.X_OK),
        }

        permission_issues = []
        if not permissions["read"]:
            permission_issues.append("read")
        if not permissions["write"]:
            permission_issues.append("write")
        if not permissions["execute"]:
            permission_issues.append("execute")

        if permission_issues:
            self.show_error(f"Permission denied: "
                            f"Unable to {' and '.join(permission_issues)} '{path}'.")
            return False

        return True

    def load_stylesheet(self):
        """
        Loads the stylesheet from an external file.

        :return: None
        """
        try:
            with open(os.path.join(self.project_root,
            "styles", "styles.qss"), "r", encoding="utf-8") as file:
                stylesheet = file.read()
                self.setStyleSheet(stylesheet)
        except (FileNotFoundError, OSError) as e:
            self.show_error(f"Error loading stylesheet: {str(e)}")

    def load_csv_on_startup(self, csv_path):
        """
        Loads an existing CSV file and populates the experiments table with its contents.
        :param csv_path: Path to the CSV file to load.
        :return: None
        """
        try:
            if os.path.exists(csv_path):
                with open(csv_path, 'r', newline='', encoding='utf-8') as file:
                    reader = csv.reader(file)
                    next(reader)
                    for row in reader:
                        if row:
                            row_position = self.experiments_table.rowCount()
                            self.experiments_table.insertRow(row_position)
                            self.experiments_table.setItem(
                                row_position, 0, QTableWidgetItem(row[0]))
                            self.experiments_table.setItem(row_position,
                            1, QTableWidgetItem(row[1]))
                self.experiments_table.resizeRowsToContents()
        except (OSError, IOError) as e:
            self.show_error(f"Error reading CSV file: {str(e)}")

    def load_yaml_on_startup(self, yaml_path):
        """
        Loads an existing YAML configuration file and updates UI elements accordingly.

        :param yaml_path: Path to the YAML file to load.
        :return: None
        """
        try:
            if os.path.exists(yaml_path):
                with open(yaml_path, 'r', encoding='utf-8') as file:
                    config = yaml.safe_load(file)
                self.basecall_checkbox.setChecked(config.get('basecall', False))
                self.dorado_model_combobox.setCurrentText(config.get('dorado_model', "fast"))
                self.unset_unsaved_changes()
        except (OSError, IOError) as e:
            self.show_error(f"Error reading YAML file: {str(e)}")
        except yaml.YAMLError as e:
            self.show_error(f"Error parsing YAML file: {str(e)}")

    def save_yaml(self):
        """
        Saves the user input to a predefined YAML file, overwriting any existing file.
        """
        try:
            basecall = self.basecall_checkbox.isChecked()
            dorado_model = self.dorado_model_combobox.currentText()
            yaml_content = (
                "# DEFAULT\n"
                "# BASECALLING\n"
                "# define parameters for a specific rule\n"
                f"basecall: {'True' if basecall else 'False'}\n"
                "# dorado model 'fast', 'hac' or 'sup'\n"
                f"dorado_model: {dorado_model}\n"
            )

            yaml_path = os.path.join(self.project_root, "config", "config.yaml")
            with open(yaml_path, 'w', encoding='utf-8') as file:
                file.write(yaml_content)

            self.log_window.append("YAML file saved successfully.")

        except (OSError, IOError) as e:
            self.show_error(f"Error saving YAML file: {str(e)}")

    def save_experiments_csv(self):
        """
        Saves the experiments data to a predefined CSV file, overwriting any existing file.
        """
        try:
            if self.experiments_table.rowCount() == 0:
                self.show_error("No experiments to save.")
                return

            csv_path = os.path.join(self.project_root, "config", "experiments.csv")
            with open(csv_path, 'w', newline='', encoding='utf-8') as file:
                writer = csv.writer(file)
                writer.writerow(["Path to Sample", "Path to Reference"])  # Write header
                for row in range(self.experiments_table.rowCount()):
                    sample_path = self.experiments_table.item(row,
                    0).text() if self.experiments_table.item(row,0) else ""
                    reference_path = self.experiments_table.item(row,
                    1).text() if self.experiments_table.item(row,1) else ""
                    writer.writerow([sample_path, reference_path])

            self.log_window.append("Experiments CSV file saved successfully.")

        except (OSError, IOError) as e:
            self.show_error(f"Error saving CSV file: {str(e)}")

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

    def set_unsaved_changes(self):
        """
        Setting unsaved changes to True
        :return:
        """
        self.unsaved_changes = True

    def unset_unsaved_changes(self):
        """
        Setting unsaved changes to False
        :return:
        """
        self.unsaved_changes = False

    def closeEvent(self, event):
        """
        Handles the window close event to manage unsaved
        changes before exiting the application.

        This method prompts the user to save any unsaved changes
        when the user attempts to close the window.
        It displays a confirmation dialog asking whether to save
        the changes, and if the user chooses to do so,
        it presents additional prompts to save the CSV and YAML
        files if applicable.

        :param event:
        :return:
        """
        # pylint: disable=C0103
        if self.unsaved_changes:
            reply = QMessageBox.question(
                self,
                'Message',
                "You have unsaved changes. Do you want to save them?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Cancel
            )
            if reply == QMessageBox.StandardButton.Yes:
                csv_reply = QMessageBox.question(
                    self,
                    'Save CSV File',
                    "Do you want to save the CSV file?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    QMessageBox.StandardButton.No
                )
                if csv_reply == QMessageBox.StandardButton.Yes:
                    self.save_experiments_csv()
                yaml_reply = QMessageBox.question(
                    self,
                    'Save YAML File',
                    "Do you want to save the YAML file?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    QMessageBox.StandardButton.No
                )
                if yaml_reply == QMessageBox.StandardButton.Yes:
                    self.save_yaml()

                self.unset_unsaved_changes()
                event.accept()
            elif reply == QMessageBox.StandardButton.No:
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()

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
        self.set_unsaved_changes()
        row_position = self.experiments_table.rowCount()
        self.experiments_table.insertRow(row_position)
        self.experiments_table.resizeRowsToContents()

    def delete_row(self):
        """
        Deletes the selected rows from the experiments table after user confirmation.
        :return: None
        """
        selected_indexes = self.experiments_table.selectedIndexes()
        selected_rows = set(index.row() for index in selected_indexes)
        if selected_rows:
            reply = QMessageBox.question(
                self,
                "Confirm Deletion",
                "Are you sure you want to delete the selected rows?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No
            )
            if reply == QMessageBox.StandardButton.Yes:
                for row in sorted(selected_rows, reverse=True):
                    self.experiments_table.removeRow(row)
                self.experiments_table.clearSelection()
                self.set_unsaved_changes()
                self.experiments_table.setCurrentCell(-1, -1)
                self.experiments_table.resizeRowsToContents()
        else:
            self.show_error("Please select a row to delete.")



if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = YamlForm()
    ex.show()
    sys.exit(app.exec())
