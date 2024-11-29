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
import signal
import subprocess
import sys
import select
import yaml
from PyQt6.QtCore import Qt  # pylint: disable=E0611
from PyQt6.QtWidgets import (QApplication, QWidget, QVBoxLayout,
                             QHBoxLayout, QLabel, QComboBox, QPushButton,
                             QCheckBox, QTableWidget, QTableWidgetItem,
                             QHeaderView, QAbstractItemView, QMessageBox,
                             QTextEdit, QSplitter, QFileDialog, QLineEdit)
from PyQt6.QtCore import QThread, pyqtSignal

class RunWorker(QThread):
    output_signal = pyqtSignal(str)
    error_signal = pyqtSignal(str)
    finished_signal = pyqtSignal()

    def __init__(self, command, parent=None,stop_requested=False):
        super().__init__(parent)
        self.command = command
        self.process = None
        self.stop_requested = stop_requested

    def run(self):
        """
        Run the shell command and capture the output in the background.
        """
        self.process = subprocess.Popen(
            ["bash", "-c", self.command],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        while True:
            rlist, _, _ = select.select([self.process.stdout, self.process.stderr], [], [], 0.1)

            if self.stop_requested:
                self._terminate_process()
                self.finished_signal.emit()
                return

            for stream in rlist:
                output = stream.readline()
                if output == '':
                    continue
                elif stream == self.process.stdout:
                    self.output_signal.emit(output.strip())
                elif stream == self.process.stderr:
                    self.error_signal.emit(output.strip())

            if self.process.poll() is not None:
                self.finished_signal.emit()
                break

    def _terminate_process(self):
        """
        Ensure the child process and all its children are terminated.
        """
        if self.process:
            self.process.terminate()
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.process.kill()

            if self.process.pid:
                try:
                    os.killpg(self.process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
            self.process = None


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

        self.worker = None
        if getattr(sys, 'frozen', False):
            self.project_root = os.path.dirname(sys.executable)
        else:
            self.project_root = os.path.abspath(os.path.join(os.path.dirname(__file__)))

        self.yaml_path = os.path.join(self.project_root, "config", "config.yaml")
        self.csv_path = os.path.join(self.project_root, "config", "experiments.csv")

        if not self.check_permissions(self.yaml_path) or not self.check_permissions(self.csv_path):
            return

        self.run_stop_button = None
        self.stop_requested = False
        self.log_window = None


        self.columns_config = []
        self.yaml_config = {}
        self.experiments_table = None

        self.left_widget = QWidget()
        self.left_widget.setObjectName("left_widget")
        self.left_layout = QVBoxLayout(self.left_widget)

        self.unsaved_changes = None

        self.widget_classes = {
            'checkbox': (QCheckBox, lambda w, v=None:
            w.setChecked(bool(v)) if v is not None else w.isChecked()),
            'combobox': (QComboBox, lambda w, v=None:
            w.setCurrentText(v) if v else w.currentText()),
            'line': (QLineEdit, lambda w, v=None:
            w.setText(v) if v else w.text()),
            'path': (QLineEdit, lambda w, v=None:
            w.setText(v) if v else w.text())
        }

        self.load_yaml_on_startup(self.yaml_path)
        self.initUI()
        self.load_csv_on_startup(self.csv_path)

        self.unset_unsaved_changes()

        self.showMaximized()
        self.running = False



    def initUI(self):
        """
        Sets up the user interface for the form.
        :return: None
        """
        # pylint: disable=C0103
        main_layout = QVBoxLayout()

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setHandleWidth(10)
        # SECTION: Left side - Controls and Log Window
        # SECTION: Controls (Run button, Basecall checkbox, Dorado model combobox, Update button)
        update_button = QPushButton("Update", self)
        self.left_layout.addWidget(update_button)
        self.create_widget_in_layout()
        self.run_stop_button = QPushButton("Run", self)
        self.run_stop_button.clicked.connect(self.toggle_run_stop)
        self.left_layout.addWidget(self.run_stop_button)

        # SECTION: Log Window
        self.log_window = QTextEdit(self)
        self.log_window.setObjectName("log_window")
        self.log_window.setReadOnly(True)
        self.left_layout.addWidget(self.log_window)

        # SECTION: Right side - Experiments Table
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)

        right_layout.addWidget(QLabel("Experiments"))
        self.experiments_table = QTableWidget(self)
        right_layout.addWidget(self.experiments_table)
        self.setup_table_columns()

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
        splitter.addWidget(self.left_widget)
        splitter.setSizes([self.width() * 2 // 3, self.width() // 3])

        main_layout.addWidget(splitter)

        self.setLayout(main_layout)
        self.setWindowTitle('Oligo-bench')
        self.load_stylesheet()

    def toggle_run_stop(self):
        """
        Toggles the Run/Stop button text and state,
         prompting to save CSV and YAML files if necessary.
        """
        if self.running:
            self.handle_stop()
        else:
            if self.experiments_table.rowCount() == 0:
                self.show_error("No cells in experiments table")
                return
            for row in range(self.experiments_table.rowCount()):
                for col in range(self.experiments_table.columnCount()):
                    item = self.experiments_table.item(row, col)
                    if item is None or not item.text().strip():
                        self.show_error("Some cells in experiments table are empty.")
                        return
            csv_reply = QMessageBox.question(
                self,
                'Save CSV File',
                "Do you want to save the CSV file before running?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Cancel
            )

            if csv_reply == QMessageBox.StandardButton.Yes:
                self.save_experiments_csv()
            if csv_reply == QMessageBox.StandardButton.Cancel:
                return
            yaml_reply = QMessageBox.question(
                self,
                'Save YAML File',
                "Do you want to save the YAML file before running?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No| QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Cancel
            )
            if yaml_reply == QMessageBox.StandardButton.Yes:
                self.save_yaml()
            if yaml_reply == QMessageBox.StandardButton.Cancel:
                return
            self.log_window.clear()
            self.run_stop_button.setText("Stop")
            self.run_stop_button.setStyleSheet("background-color: #9b3438; color: white;")
            self.unset_unsaved_changes()
            self.running = True
            self.stop_requested = False
            self.start_run()

    def start_run(self):
        """
        Runs the shell script (run.sh) and captures the output.
        Handles errors and displays localized, specific error messages.
        """
        self.log_window.clear()
        self.stop_requested = False

        try:
            path = os.path.join(self.project_root, 'run.sh')

            try:
                conda_bin_path = subprocess.check_output(["which", "conda"]).decode().strip()
            except subprocess.CalledProcessError:
                raise FileNotFoundError("Conda is not installed or not found in the system path. Please install Conda.")

            conda_dir = os.path.dirname(conda_bin_path)
            conda_init_path = os.path.join(conda_dir, "..", "etc", "profile.d", "conda.sh")

            if not os.path.exists(conda_init_path):
                raise FileNotFoundError(
                    f"Conda initialization script not found at {conda_init_path}. Please check your Conda installation.")

            env_path = 'snakemake'
            command = (
                f"source {conda_init_path} && "
                f"conda activate {env_path} && "
                f"bash {path}"
            )

            self.worker = RunWorker(command)
            self.worker.output_signal.connect(self.handle_output)
            self.worker.error_signal.connect(self.handle_error)
            self.worker.finished_signal.connect(self.finish_run)
            self.worker.start()

        except FileNotFoundError as e:
            self.show_error(
                f"Error: {str(e)}\nPlease ensure all necessary files and directories exist and are properly configured.")
        except subprocess.CalledProcessError as e:
            self.show_error(
                f"Error: A command execution failed. Please check the command and ensure Conda is properly initialized.\nError details: {e.output.decode()}")

    def handle_output(self, output):
        """
        Handles the output from the worker thread and updates the log window.
        """
        self.log_window.append(f"<font color='#b4b2af'><pre>{output}</pre></font>")

    def handle_error(self, error_output):
        """
        Handles the error output from the worker thread and updates the log window.
        """
        self.log_window.append(f"<font color='#9b3438'><pre>Error: {error_output}</pre></font>")

    def handle_stop(self):
        """
        Handles stopping the execution and updates the UI.
        """
        self.stop_requested = True
        self.run_stop_button.setText("Stopping...")
        self.run_stop_button.setStyleSheet("background-color: yellow; color: black;")
        if self.worker:
            self.worker.stop_requested = True
        self.log_window.append("<font color='#9b3438'>Execution was stopped.</font>")
        self.run_stop_button.setText("Run")
        self.run_stop_button.setStyleSheet("")


    def finish_run(self):
        """
        Resets the state after the script finishes.
        """
        self.worker = None
        self.running = False
        self.run_stop_button.setText("Run")
        self.run_stop_button.setStyleSheet("")
        self.log_window.append("<font color='yellow'>Run finished.</font>")

    def handle_cell_double_click(self, row, column):
        """
        Called when a cell is double-clicked.
        :param row: Row of the cell
        :param column: Column of the cell
        """
        column_type = self.columns_config[column].get("type")
        if column_type == "plaintxt":
            item = self.experiments_table.item(row, column)
            if not item:
                item = QTableWidgetItem()
                self.experiments_table.setItem(row, column, item)
                self.set_unsaved_changes()
            self.experiments_table.editItem(item)
            self.set_unsaved_changes()
        else:
            self.table_cell_open_file(row,column)

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
            if not os.path.exists(csv_path):
                self.show_error("CSV path does not exist.")
                return
            with open(csv_path, 'r', newline='', encoding='utf-8') as file:
                reader = csv.reader(file)
                next(reader, None)
                for row in reader:
                    if not row:
                        continue
                    row_position = self.experiments_table.rowCount()
                    self.experiments_table.insertRow(row_position)
                    for col_index, col_config in enumerate(self.columns_config):
                        cell_value = row[col_index] if col_index < len(row) else ""
                        self.experiments_table.setItem(row_position, col_index, QTableWidgetItem(cell_value))
        except (OSError, IOError) as e:
            self.show_error(f"Error reading CSV file: {str(e)}")
            return
        self.experiments_table.resizeRowsToContents()

    def load_yaml_on_startup(self, yaml_path):
        """
        Loads an existing YAML file.
        :param yaml_path:
        :return:
        """
        try:
            if os.path.exists(yaml_path):
                with open(yaml_path, 'r', encoding='utf-8') as file:
                    config = yaml.safe_load(file)
                self.yaml_config = config.get('basecalling', {})
                self.columns_config = config.get('experiments', [])
                self.unset_unsaved_changes()
        except (OSError, IOError) as e:
            self.show_error(f"Error reading YAML file: {str(e)}")
        except yaml.YAMLError as e:
            self.show_error(f"Error parsing YAML file: {str(e)}")

    def setup_table_columns(self):
        """
        Sets up columns in the experiments table based on the YAML configuration.
        """
        if self.columns_config:
            self.experiments_table.setColumnCount(len(self.columns_config))
            self.experiments_table.setHorizontalHeaderLabels([col["column"] for col in self.columns_config])
            self.experiments_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
            self.experiments_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
            header = self.experiments_table.horizontalHeader()
            header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
            for col_index in range(len(self.columns_config)):
                header.resizeSection(col_index, 300)
            self.experiments_table.cellDoubleClicked.connect(self.handle_cell_double_click)

    def save_yaml(self):
        """
        Saves only the updated 'value' fields back to the YAML file, preserving the rest of the configuration.
        :return:
        """
        try:
            with open(self.yaml_path, 'r', encoding='utf-8') as file:
                current_config = yaml.safe_load(file)
        except (OSError, IOError) as e:
            self.show_error(f"Error opening YAML file: {str(e)}")
            return
        except yaml.YAMLError as e:
            self.show_error(f"Error reading YAML data: {str(e)}")
            return

        if 'basecalling' not in current_config:
            current_config['basecalling'] = {}

        for widget_name, widget_data in self.yaml_config.items():
            widget_type = widget_data.get('type')
            object_name = f"{widget_name}_{widget_type}"
            if widget_type in self.widget_classes and widget_name in current_config['basecalling']:
                widget_class, value_getter = self.widget_classes[widget_type]
                try:
                    widget = self.findChild(widget_class, object_name)
                    if widget:
                        current_config['basecalling'][widget_name]['value'] = value_getter(widget)
                except AttributeError as e:
                    self.show_error(f"Error accessing widget {object_name}: {str(e)}")

        try:
            with open(self.yaml_path, 'w', encoding='utf-8') as file:
                yaml.dump(current_config, file, default_flow_style=False)
        except (OSError, IOError) as e:
            self.show_error(f"Error writing to YAML file: {str(e)}")
            return
        except yaml.YAMLError as e:
            self.show_error(f"Error handling YAML data on save: {str(e)}")
            return

        self.log_window.append("YAML file saved successfully.")
        self.unset_unsaved_changes()

    def save_experiments_csv(self):
        """
        Saves the experiments data to a predefined CSV file, overwriting any existing file.
        """

        if self.experiments_table.rowCount() == 0:
            self.show_error("No experiments to save.")
            return
        try:
            with open(self.csv_path, 'w', newline='', encoding='utf-8') as file:
                writer = csv.writer(file)
                headers = [col["column"] for col in self.columns_config]
                writer.writerow(headers)
                for row in range(self.experiments_table.rowCount()):
                    row_data = [
                        self.experiments_table.item(row, col).text() if self.experiments_table.item(row, col) else ""
                        for col in range(self.experiments_table.columnCount())
                    ]
                    writer.writerow(row_data)
        except (OSError, IOError) as e:
            self.show_error(f"Error saving CSV file: {str(e)}")
        self.log_window.append("Experiments CSV file saved successfully.")

    def table_cell_open_file(self, row, column):
        """
        Opens a file dialog to select a file or directory based on the column type
        specified in the YAML configuration and updates the specified cell.

        :param row: The row index of the cell that was double-clicked
        :param column: The column index of the cell that was double-clicked
        :return: None
        """
        if not self.columns_config or column >= len(self.columns_config):
            return
        column_type = self.columns_config[column].get("type")
        cell_item = self.experiments_table.item(row, column)
        current_text = cell_item.text() if cell_item else ""
        selected_path = self.open_file_dialog(column_type, current_text)
        if not cell_item:
            cell_item = QTableWidgetItem()
            self.experiments_table.setItem(row, column, cell_item)
        cell_item.setText(selected_path)

    def open_file_dialog(self, file_type: str, old_value: str):
        """
        Opens a file dialog to select a file or directory based on the type of file.
        If the selected path is empty, it returns the old value instead.

        :param file_type: str
        :param old_value: str (the previous path to fall back to if nothing is selected)
        :return: selected_path: str
        """
        if file_type == "folder":
            selected_path = QFileDialog.getExistingDirectory(self, "Select Folder", "")
        else:
            file_filter = f"{file_type.upper()} Files (*.{file_type.lower()})"
            selected_path, _ = QFileDialog.getOpenFileName(self, f"Select {file_type.upper()} File", "", file_filter)

        if selected_path:
            self.set_unsaved_changes()
            return selected_path
        else:
            return old_value

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

    def create_widget_in_layout(self):
        """
        Creates a widget dynamically based on configuration and adds it to the provided layout.
        :return: None
        """
        for widget_name, widget_config in self.yaml_config.items():
            if widget_name and isinstance(widget_config, dict):
                widget_type = widget_config.get('type')
                widget_label = widget_config.get('label')
                widget_value = widget_config.get('value', None)
                object_name = f"{widget_name}_{widget_type}"

                widget = None
                if widget_type in self.widget_classes:
                    widget_class, set_value = self.widget_classes[widget_type]
                    try:
                        widget = widget_class()
                        widget.setObjectName(object_name)
                        set_value(widget, widget_value)
                        if widget_type == 'checkbox':
                            widget.stateChanged.connect(self.set_unsaved_changes)
                        elif widget_type == 'combobox':
                            widget_options = widget_config.get('options', [])
                            if widget_options:
                                widget.addItems(widget_options)
                            widget.currentIndexChanged.connect(self.set_unsaved_changes)
                        elif widget_type == 'line':
                            widget.textChanged.connect(self.set_unsaved_changes)
                        elif widget_type == 'path':
                            widget_button: QPushButton = QPushButton("Open File", self)
                            widget_button.setObjectName(f"{widget_name}_{widget_type}_button")
                            widget_button.clicked.connect(
                                lambda: widget.setText(
                                    self.open_file_dialog(widget_config.get('file_type'), widget.text())))
                            widget.textChanged.connect(self.set_unsaved_changes)

                    except Exception as e:
                        self.show_error(f"Error creating widget {object_name}: {str(e)}")

                if widget:
                    try:
                        name_label = QLabel(f"{widget_label}:")
                        temp_widget = QWidget()
                        widget_layout = QHBoxLayout(temp_widget)
                        widget_layout.addWidget(name_label)
                        widget_layout.addWidget(widget)
                        temp_widget.setObjectName(f"{object_name}_layout")
                        name_label.setObjectName(f"{object_name}_name_label")
                        if widget_type == 'path':
                            widget_layout.addWidget(widget_button)
                        self.left_layout.addWidget(temp_widget)
                    except Exception as e:
                        self.show_error(f"Error adding widget {object_name} to layout: {str(e)}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = YamlForm()
    ex.show()
    sys.exit(app.exec())
