"""
YamlForm Application

This module defines a PyQt6 application that provides
a graphical user interface (GUI) for managing YAML and CSV files
related to experiments. The application allows users
to configure base calling options, select dorado models, and
manage experiment paths through a table interface.
"""
import csv
import logging
import os
import shutil
import signal
import subprocess
import sys


import PyQt6.QtCore
import requests
import yaml
from PyQt6.QtCore import Qt  # pylint: disable=E0611
from PyQt6.QtWidgets import (QApplication, QWidget, QVBoxLayout,
                             QHBoxLayout, QLabel, QComboBox, QPushButton,
                             QCheckBox, QTableWidget, QTableWidgetItem,
                             QHeaderView, QAbstractItemView, QMessageBox,
                             QTextEdit, QSplitter, QFileDialog, QLineEdit, QSpinBox)
from PyQt6.QtCore import QThread, pyqtSignal

log_dir = os.path.join(os.path.dirname(__file__), "logs")
os.makedirs(log_dir, exist_ok=True)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
log_file_path = os.path.join(log_dir, "execution.log")
handler = logging.FileHandler(log_file_path, mode='w')
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


class UpdateWorker(QThread):
    output_signal = pyqtSignal(str)
    error_signal = pyqtSignal(str)
    finished_signal = pyqtSignal()
    def __init__(self,project_root,backup_path, parent=None):
        super().__init__(parent)
        self.project_root= project_root
        self.backup_path = backup_path


    def run(self):
        """
        Updates the project directory with the latest contents from a remote repository.

        This function performs the following steps:

        1. Moves the `config` and oligo environment folders to temporary locations outside the project directory.
        2. Deletes all files and folders in the project directory except for the `logs` folder.
        3. Downloads a zip archive from a specified URL and extracts its contents.
        4. Moves the contents of the extracted `oligo-bench-main` folder to the project directory.
        5. Restores the `oligo` folder from the temporary location back to the project directory.
        6. Replaces the `experiments.csv` file in the `config` folder with the version from the temporary location.
        7. Transfers configuration values from a temporary `config.yaml` file to the restored `config.yaml` file.
        8. Cleans up temporary folders and ensures the project directory is updated.
        :return:
        """

        zip_url = "https://github.com/genomika-lt/oligo-bench/archive/refs/heads/main.zip"
        zip_path = os.path.join(self.project_root, "main.zip")

        oligo_path = os.path.join(self.project_root, "oligo")
        config_folder_path = os.path.join(self.project_root, "config")
        dorado_path = os.path.join(self.project_root, "dorado")
        temp_oligo_path = os.path.join(self.backup_path, "oligo")
        temp_config_path = os.path.join(self.backup_path, "config")
        temp_dorado_path = os.path.join(self.backup_path, "dorado")

        self.output_signal.emit("Starting download process...")
        logger.info("Starting download process...")

        if os.path.exists(self.backup_path):
            shutil.rmtree(self.backup_path)

        os.makedirs(self.backup_path, exist_ok=True)
        for item in os.listdir(self.project_root):
            item_path = os.path.join(self.project_root, item)
            if os.path.basename(item_path) != "logs":
                shutil.move(item_path, self.backup_path)
        self.output_signal.emit(f"Moved project contents to rollback backup at {self.backup_path}")
        logger.info(f"Moved project contents to rollback backup at {self.backup_path}")

        try:
            self.output_signal.emit(f"Downloading archive from {zip_url}")
            logger.info(f"Downloading archive from {zip_url}")
            response = requests.get(zip_url, stream=True)
            response.raise_for_status()
            with open(zip_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
        except requests.RequestException as e:
            self.error_signal.emit(f"Failed to download the zip archive: {e}")
            logger.error(f"Failed to download the zip archive: {e}")
            raise RuntimeError(f"Failed to download the zip archive: {e}")
        self.output_signal.emit(f"Downloaded archive to {zip_path}")
        logger.info(f"Downloaded archive to {zip_path}")

        try:
            os.system("unzip main.zip")
        except Exception as e:
            self.error_signal.emit(f"Error unpacking zip file: {e}")
            logger.error(f"Error unpacking zip file: {e}")
            raise e
        self.output_signal.emit("Extracted archive contents")
        logger.info("Extracted archive contents")

        os.remove(zip_path)
        self.output_signal.emit(f"Removed zip file in {zip_path}")
        logger.info(f"Removed zip file in {zip_path}")

        extracted_folder = os.path.join(self.project_root, "oligo-bench-main")
        for item in os.listdir(extracted_folder):
            shutil.move(os.path.join(extracted_folder, item), self.project_root)
        shutil.rmtree(extracted_folder)
        self.output_signal.emit("Moved contents from downloaded zip to project root")
        logger.info("Moved contents from downloaded zip to project root")

        shutil.move(temp_oligo_path, oligo_path)
        self.output_signal.emit("Moved oligo environment folder from temporary path to project root")
        logger.info("Moved oligo environment folder from temporary path to project root")

        shutil.move(temp_dorado_path, dorado_path)
        self.output_signal.emit("Moved dorado folder from temporary path to project root")
        logger.info("Moved dorado folder from temporary path to project root")

        experiments_csv_path = os.path.join(config_folder_path, "experiments.csv")
        os.remove(experiments_csv_path)
        self.output_signal.emit(f"Removed experiments default csv file in {experiments_csv_path}")
        logger.info(f"Removed experiments default csv file in {experiments_csv_path}")

        shutil.copy(os.path.join(temp_config_path, "experiments.csv"), experiments_csv_path)
        self.output_signal.emit("Replaced default experiments.csv file with saved data")
        logger.info("Replaced default experiments.csv file with saved data")

        temp_config_yaml = os.path.join(temp_config_path, "config.yaml")
        self.transfer_values(temp_config_yaml, os.path.join(config_folder_path, "config.yaml"))

        shutil.rmtree(self.backup_path)
        self.output_signal.emit("Removed backup folder for rolling back changes")
        logger.info("Removed backup folder for rolling back changes")

        self.output_signal.emit("Project update completed successfully")
        logger.info("Project update completed successfully")
        self.finished_signal.emit()
        
        
    def transfer_values(self,source_file: str, target_file: str):
        """
        Transfers the `value` field of every widget in the `parameters` section
        from the source YAML file to the target YAML file.

        Args:
            source_file (str): Path to the source YAML file.
            target_file (str): Path to the target YAML file.
        """
        try:
            with open(source_file, 'r') as src:
                source_data = yaml.safe_load(src)
        except (FileNotFoundError, IOError, yaml.YAMLError) as e:
            self.error_signal.emit(f"Error reading source file '{source_file}': {e}")
            logger.error(f"Error reading source file '{source_file}': {e}")
            raise

        try:
            with open(target_file, 'r') as tgt:
                target_data = yaml.safe_load(tgt)
        except (FileNotFoundError, IOError, yaml.YAMLError) as e:
            self.error_signal.emit(f"Error reading target file '{target_file}': {e}")
            logger.error(f"Error reading target file '{target_file}': {e}")
            raise

        if 'parameters' not in source_data or 'parameters' not in target_data:
            self.error_signal.emit("Validation error: 'parameters' section missing")
            logger.error("Validation error: 'parameters' section missing")
            raise KeyError("Both source and target files must contain a 'parameters' section.")

        try:
            for key, widget in source_data['parameters'].items():
                if 'value' in widget and key in target_data['parameters']:
                    target_data['parameters'][key]['value'] = widget['value']
        except (AttributeError, TypeError) as e:
            self.error_signal.emit(f"Error processing 'parameters' sections: {e}")
            logger.error(f"Error processing 'parameters' sections: {e}")
            raise

        try:
            with open(target_file, 'w') as tgt:
                yaml.safe_dump(target_data, tgt)
        except (IOError, yaml.YAMLError) as e:
            self.error_signal.emit(f"Error writing to target file '{target_file}': {e}")
            logger.error(f"Error writing to target file '{target_file}': {e}")
            raise

        self.output_signal.emit(f"Successfully transferred values from '{source_file}' to '{target_file}'.")
        logger.info(f"Successfully transferred values from '{source_file}' to '{target_file}'.")

            
class RunWorker(QThread):
    output_signal = pyqtSignal(str)
    error_signal = pyqtSignal(str)
    finished_signal = pyqtSignal()
    def __init__(self, command,project_root, parent=None, stop_requested=False):
        super().__init__(parent)
        self.command = command
        self.process = None
        self.stop_requested = stop_requested
        self.project_root= project_root


    def run(self):
        """
        Run the shell command and capture the output in the background.
        """
        logger.info(f"Starting command: {self.command}")

        self.process = subprocess.Popen(
            ["bash", "-c", self.command],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        while True:
            if self.stop_requested:
                logger.info("Stop requested. Terminating process.")
                self._terminate_process()
                self.finished_signal.emit()
                return
            if self.process.poll() is not None:
                logger.info(f"Process finished with return code {self.process.returncode}")
                self.finished_signal.emit()
                break
            stderr_output = self.process.stderr.readline()
            if self._is_real_error(stderr_output.strip()):
                self.error_signal.emit(stderr_output.strip())
                logger.error(stderr_output.strip())
            else:
                self.output_signal.emit(stderr_output.strip())
                logger.info(stderr_output.strip())


    def _is_real_error(self, line):
        """
        Determine if a line from stderr represents a real error.
        """
        excluded_keywords = ["reason:", "input:", "log:","resources:","output:"]
        line_lower = line.lower()
        if "error:" in line_lower:
            if not any(excluded in line_lower for excluded in excluded_keywords):
                return True
        return False


    def _terminate_process(self):
        """
        Ensure the child process and all its children are terminated.
        """
        if self.process:
            logger.info("Terminating process.")
            self.process.terminate()
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                logger.error("Process timeout. Force killing the process.")
                self.process.kill()

            if self.process.pid:
                try:
                    os.killpg(self.process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    logger.error(f"Failed to terminate process group: {self.process.pid}")
            self.process = None


def load_stylesheet(widget):
    """
    Loads the stylesheet from an external file.

    :return: None
    """
    if getattr(sys, 'frozen', False):
        css_path = os.path.join(sys._MEIPASS, "styles", "styles.qss")
    else:
        css_path = os.path.join(os.path.dirname(__file__), "styles",
                                "styles.qss")
    try:
        with open(css_path, "r", encoding="utf-8") as file:
            stylesheet = file.read()
            widget.setStyleSheet(stylesheet)
    except (FileNotFoundError, OSError) as e:
        show_error(f"Error loading stylesheet: {str(e)}")


def show_error(widget, message):
    """
    Displays an error message in a popup window.
    :param widget: QWidget(parent) of QMessageBox.
    :param message: message to display
    :return: None
    """
    QMessageBox.critical(widget, "Error", message)


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
        if getattr(sys, 'frozen', False):
            self.project_root = os.path.dirname(sys.executable)
        else:
            self.project_root = os.path.abspath(os.path.join(os.path.dirname(__file__)))
        self.yaml_path = os.path.join(self.project_root, "config", "config.yaml")
        self.csv_path = os.path.join(self.project_root, "config", "experiments.csv")

        self.run_stop_button = None
        self.update_button = None
        self.delete_row_button = None
        self.add_row_button = None
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
            w.setText(v) if v else w.text()),
            'integer': (QSpinBox, lambda w, v=None:
            w.setValue(v) if v else w.value()),
        }
        perm_yaml_path=self.check_permissions(self.yaml_path)
        perm_csv_path = self.check_permissions(self.csv_path)
        if not perm_yaml_path:
            logger.error(f"{self.yaml_path} does not have required permissions.")
            show_error(self,f"{self.yaml_path} does not have required permissions.")
        if not perm_csv_path:
            logger.error(f"{self.csv_path} does not have required permissions.")
            show_error(self,f"{self.csv_path} does not have required permissions.")
        if perm_yaml_path and perm_csv_path:
            self.load_yaml_on_startup(self.yaml_path)
            self.initUI()
            self.load_csv_on_startup(self.csv_path)

        self.worker = None
        self.update_worker = None
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
        # SECTION: Controls (Run button, Basecall checkbox, Dorado model combobox, Settings button)
        self.settings_button = QPushButton("Settings", self)
        self.settings_button.clicked.connect(self.open_settings_window)
        self.left_layout.addWidget(self.settings_button)
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
        self.experiments_table.setFocusPolicy(PyQt6.QtCore.Qt.FocusPolicy.NoFocus)
        right_layout.addWidget(self.experiments_table)
        self.setup_table_columns()

        # Buttons for adding and deleting rows
        table_buttons_layout = QHBoxLayout()

        self.add_row_button = QPushButton("Add", self)
        self.add_row_button.clicked.connect(self.add_row)
        table_buttons_layout.addWidget(self.add_row_button)

        self.delete_row_button = QPushButton("Delete", self)
        self.delete_row_button.clicked.connect(self.delete_row)
        table_buttons_layout.addWidget(self.delete_row_button)

        right_layout.addLayout(table_buttons_layout)

        splitter.addWidget(right_widget)
        splitter.addWidget(self.left_widget)
        splitter.setSizes([self.width() * 2 // 3, self.width() // 3])

        main_layout.addWidget(splitter)

        self.setLayout(main_layout)
        self.setWindowTitle('Oligo-bench')
        load_stylesheet(self)


    def toggle_run_stop(self):
        """
        Toggles the Run/Stop button text and state,
         prompting to save CSV and YAML files if necessary.
        """
        if self.running:
            self.handle_stop()
            return
        if self.experiments_table.rowCount() == 0:
            show_error(self,"No cells in experiments table")
            return
        for row in range(self.experiments_table.rowCount()):
            for col in range(self.experiments_table.columnCount()):
                item = self.experiments_table.item(row, col)
                if item is None or not item.text().strip():
                    show_error(self,"Some cells in experiments table are empty.")
                    return
        if self.unsaved_changes:
            reply = QMessageBox.question(
                self,
                'Save CSV and YAML files',
                "Do you want to save the CSV and YAML files before running?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Cancel
            )
            if reply == QMessageBox.StandardButton.Yes:
                self.save_experiments_csv()
                self.save_yaml()
            if reply == QMessageBox.StandardButton.Cancel:
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

        path = os.path.join(self.project_root, './run.sh')
        if not path:
            show_error(self,f"File {path} does not exist")
            logger.error(f"File {path} does not exist")
            return

        try:
            subprocess.check_output(["which", "conda"]).decode().strip()
        except subprocess.CalledProcessError:
            raise FileNotFoundError("Conda is not installed or not found in the system path. Please install Conda.")

        env_path = os.path.join(self.project_root,'oligo')
        command = f"conda run --live-stream -p {env_path} bash {path}"

        self.worker = RunWorker(command,self.project_root)
        self.worker.output_signal.connect(self.handle_output)
        self.worker.error_signal.connect(self.handle_error)
        self.worker.finished_signal.connect(self.finish_run)
        self.worker.start()


    def open_settings_window(self):
        self.settings_window = SettingsWindow(self.project_root)

        self.run_stop_button.setEnabled(False)
        self.add_row_button.setEnabled(False)
        self.delete_row_button.setEnabled(False)

        if self.running:
            self.handle_stop()

        self.settings_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self.settings_window.destroyed.connect(self.on_settings_window_closed)
        self.settings_window.show()

    def on_settings_window_closed(self):
        self.run_stop_button.setEnabled(True)
        self.add_row_button.setEnabled(True)
        self.delete_row_button.setEnabled(True)


    def handle_output(self, output):
        """
        Handles the output from the worker thread and updates the log window.
        """
        self.log_window.append(f"<font color='#49A078'><pre>{output}</pre></font>")


    def handle_error(self, error_output):
        """
        Handles the error output from the worker thread and updates the log window.
        """
        self.log_window.append(f"<font color='#9b3438'><pre>Error: {error_output}</pre></font>")
        logger.info(error_output)


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
            show_error(self,f"File not found: {path}")
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
            show_error(self,f"Permission denied: "
                            f"Unable to {' and '.join(permission_issues)} '{path}'.")
            return False

        return True


    def load_csv_on_startup(self, csv_path):
        """
        Loads an existing CSV file and populates the experiments table with its contents.
        :param csv_path: Path to the CSV file to load.
        :return: None
        """
        try:
            if not os.path.exists(csv_path):
                show_error(self,"CSV path does not exist.")
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
            show_error(self,f"Error reading CSV file: {str(e)}")
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
                self.yaml_config = config.get('parameters', {})
                self.columns_config = config.get('experiments', [])
                self.unset_unsaved_changes()
        except (OSError, IOError) as e:
            show_error(self,f"Error reading YAML file: {str(e)}")
        except yaml.YAMLError as e:
            show_error(self,f"Error parsing YAML file: {str(e)}")

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
            show_error(self,f"Error opening YAML file: {str(e)}")
            return
        except yaml.YAMLError as e:
            show_error(self,f"Error reading YAML data: {str(e)}")
            return

        if 'parameters' not in current_config:
            current_config['parameters'] = {}

        for widget_name, widget_data in self.yaml_config.items():
            widget_type = widget_data.get('type')
            object_name = f"{widget_name}_{widget_type}"
            if widget_type in self.widget_classes and widget_name in current_config['parameters']:
                widget_class, value_getter = self.widget_classes[widget_type]
                try:
                    widget = self.findChild(widget_class, object_name)
                    if widget:
                        current_config['parameters'][widget_name]['value'] = value_getter(widget)
                except AttributeError as e:
                    show_error(self,f"Error accessing widget {object_name}: {str(e)}")

        try:
            with open(self.yaml_path, 'w', encoding='utf-8') as file:
                yaml.dump(current_config, file, default_flow_style=False)
        except (OSError, IOError) as e:
            show_error(self,f"Error writing to YAML file: {str(e)}")
            return
        except yaml.YAMLError as e:
            show_error(self,f"Error handling YAML data on save: {str(e)}")
            return

        self.log_window.append("YAML file saved successfully.")
        self.unset_unsaved_changes()

    def save_experiments_csv(self):
        """
        Saves the experiments data to a predefined CSV file, overwriting any existing file.
        """

        if self.experiments_table.rowCount() == 0:
            show_error(self,"No experiments to save.")
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
            show_error(self,f"Error saving CSV file: {str(e)}")
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

        if not selected_path:
            return old_value

        self.set_unsaved_changes()
        return selected_path


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
                self.save_experiments_csv()
                self.save_yaml()
                self.unset_unsaved_changes()
                event.accept()
            elif reply == QMessageBox.StandardButton.No:
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()


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
        if not selected_rows:
            show_error(self,"Please select a row to delete.")
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
                            widget.setCurrentText(widget_value)
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
                        elif widget_type == 'integer':
                            widget.setMinimum(widget_config.get('min_value'))
                            widget.setMaximum(widget_config.get('max_value'))
                            widget.valueChanged.connect(self.set_unsaved_changes)

                    except Exception as e:
                        show_error(self,f"Error creating widget {object_name}: {str(e)}")

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
                        temp_widget.setStyleSheet("  background-color: #302c34;color:#b4b2af;")
                        widget.setStyleSheet("  background-color: white; color:#333333;")
                        self.left_layout.addWidget(temp_widget)
                    except Exception as e:
                        show_error(self,f"Error adding widget {object_name} to layout: {str(e)}")


class SettingsWindow(QWidget):
    def __init__(self, project_root):
        super().__init__()
        self.setWindowTitle("Settings")
        self.setGeometry(400, 400, 500, 300)
        self.project_root = project_root
        self.original_save_path = self.get_current_save_path()
        self.update_worker = None
        load_stylesheet(self)

        layout = QVBoxLayout(self)

        self.version_label = QLabel("Version: 1.2")

        self.update_button = QPushButton("Update", self)
        self.update_button.clicked.connect(self.update_qc)

        update_widget = QWidget()
        update_layout = QHBoxLayout(update_widget)
        update_layout.addWidget(self.version_label)
        update_layout.addWidget(self.update_button)

        layout.addWidget(update_widget)

        self.save_path_label = QLabel("Save Path:")
        self.save_path_input = QLineEdit(self.original_save_path, self)
        self.save_path_input.setReadOnly(True)
        self.save_path_input.textChanged.connect(self.check_changes)
        self.path_button = QPushButton("Choose path", self)
        self.path_button.clicked.connect(self.open_dir)

        self.apply_button = QPushButton("Apply Changes", self)
        self.apply_button.setEnabled(False)
        self.apply_button.clicked.connect(self.apply_changes)

        temp_widget = QWidget()
        path_layout = QHBoxLayout(temp_widget)
        path_layout.addWidget(self.save_path_label)
        path_layout.addWidget(self.save_path_input)
        path_layout.addWidget(self.path_button)

        layout.addWidget(temp_widget)

        layout.addWidget(self.apply_button)
        self.setLayout(layout)

    def get_current_save_path(self):
        """Retrieve the current save_path from the YAML file."""
        yaml_path = os.path.join(self.project_root, "config", "config.yaml")
        try:
            with open(yaml_path, "r") as f:
                settings = yaml.safe_load(f)
                return settings["settings"]["save_path"]
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load save path: {e}")

    def check_changes(self):
        """Enable or disable the Apply button based on changes to the save path."""
        current_text = self.save_path_input.text()
        self.apply_button.setEnabled(current_text != self.original_save_path)

    def apply_changes(self):
        yaml_path = os.path.join(self.project_root, "config", "config.yaml")
        try:
            with open(yaml_path, "r") as f:
                settings = yaml.safe_load(f)

            settings["settings"]["save_path"] = self.save_path_input.text()

            with open(yaml_path, "w") as f:
                yaml.safe_dump(settings, f)

            QMessageBox.information(self, "Success", "Settings have been updated.")
            self.original_save_path = self.save_path_input.text()
            self.apply_button.setEnabled(False)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to update settings: {e}")


    def update_qc(self):
        reply = QMessageBox.warning(
            self,
            'Update Project',
            "Do you really to update project? Any execution process during installation will be stopped",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
        )
        if reply == QMessageBox.StandardButton.Yes:
            self.update_button.setEnabled(False)
            self.apply_button.setEnabled(False)
            backup_path = os.path.join("/tmp/oligo_bench_temp")
            try:
                self.update_worker = UpdateWorker(self.project_root, backup_path)
                self.update_worker.error_signal.connect(self.handle_error)
                self.update_worker.finished_signal.connect(self.finish_update)
                self.update_worker.start()
            except Exception as e:
                self.handle_error(f"Error during project update: {e}. Rolling back changes.")
                for item in os.listdir(backup_path):
                    shutil.move(os.path.join(backup_path, item), self.project_root)
                shutil.rmtree(backup_path)
                self.handle_output("Rollback completed. Project restored to its original state.")
                raise e
        self.update_button.setEnabled(True)
        self.apply_button.setEnabled(True)

    def open_dir(self):
        """
        Opens a file dialog to select a directory and sets result to save_path in YAML file
        :return: None
        """
        current_text = self.save_path_input.text()
        selected_path = self.open_file_dialog(current_text)
        self.save_path_input.setText(selected_path)

    def open_file_dialog(self, old_value: str):
        """
        Opens a file dialog to select a directory.
        If the selected path is empty, it returns the old value instead.

        :param old_value: str (the previous path to fall back to if nothing is selected)
        :return: selected_path: str
        """
        selected_path = QFileDialog.getExistingDirectory(self, "Select Folder", "")
        if not selected_path:
            return old_value
        return selected_path

    def handle_error(self, message):
        QMessageBox.critical(self, "Error", message)

    def finish_update(self):
        QMessageBox.information(self, "Update", "Project update completed successfully.")
        self.update_worker = None
        self.update_button.setEnabled(True)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = YamlForm()
    ex.show()
    sys.exit(app.exec())
